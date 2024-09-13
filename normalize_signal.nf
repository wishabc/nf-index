params.sample_weights = ""

def non_required_arg(value, key) {
    return value ? "${key} ${value}": ""
}



process extract_from_anndata {
    conda params.conda

    input:
        path anndata
    
    output:
        tuple path("binary.matrix.npy"), path("counts.matrix.npy"), path(samples_order), path(masterlist)

    script:
    samples_order = "samples_order.txt"
    masterlist = "masterlist.no_header.bed"
    """
    python3 $moduleDir/bin/convert_to_anndata/extract_from_anndata.py \
        ${anndata} \
        ${masterlist} \
        ${samples_order} \
        --extra_layers binary,counts \
        --dhs_mask_name final_qc_passing_dhs
    """
}
process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}"

	input:
        path peaks_matrix
		path signal_matrix
		path norm_params, stageAs: "params/*"

	output:
		path "${prefix}.scale_factors.mean_normalized.npy", emit: scale_factors
        path "${prefix}.log_differences.npy", emit: log_diffs
		tuple path("${prefix}.lowess_params.npz"), path("${prefix}.lowess_params.json"), emit: model_params
        path "${prefix}.*.pdf", emit: normalization_qc
		
	script:
    pref = 'normalized.only_autosomes.filtered'
    save_dir = 'normalization'
	prefix = "${save_dir}/${pref}"
	n = norm_params.size() == 2 ? file(norm_params[0]) : ""
	normalization_params = n ? "--model_params params/${n.baseName}" : ""
	"""
	normalize-matrix ${signal_matrix} \
        ${save_dir} \
        --prefix ${pref} \
		--jobs ${task.cpus} \
        --normalization_type quantile_lowess \
		${non_required_arg(params.sample_weights, '--weights')} \
		${normalization_params}

    qc-normalization ${save_dir} \
        --prefix ${pref} \
        --n_samples 10 \
        --peak_matrix ${peaks_matrix}
	"""
}

process deseq2 {
	conda params.conda
	publishDir "${params.outdir}/normalization", pattern: "${prefix}*.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}*.RDS"
	label "bigmem"

	input:
		tuple path(scale_factors), path(signal_matrix), path(samples_order)
		path norm_params, stageAs: "params/*"

	output:
		path "${prefix}*.npy", emit: matrix
		path "${prefix}*.RDS", emit: model_params

	script:
	prefix = "deseq_normalized.only_autosomes.filtered"
	normalization_params = norm_params.name != "params/empty.params" ? norm_params : ""
	"""
	Rscript $moduleDir/bin/deseq2.R \
		${signal_matrix} \
		${scale_factors} \
		${samples_order} \
		${params.samples_file} \
		${prefix} \
		${normalization_params}
	"""
}


workflow normalizeMatrix {
	take:
		matrices  // binary_matrix, count_matrix, samples_order, masterlist
		normalization_params
	main:

		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ }
			| ifEmpty(file("empty.params"))
			| collect(sort: true)
        
        deseq_params = normalization_params
			| filter { it.name =~ /params\.RDS/ }
			| ifEmpty(file("empty.params"))

		sf = normalize_matrix(
            matrices.map(it -> tuple(it[0], it[1])),
            lowess_params
        ).scale_factors
            | combine(matrices.map(it -> tuple(it[1], it[2])))

		out = deseq2(sf, deseq_params).matrix

        //h5_file = convert_to_h5(binary_matrix, out, samples_order)

	emit:
		out
}

// Re-use existing model from other run (e.g. different samples)
workflow existingModel {
    params.template_run_dir = "/path/to/previous/run"
    if (!file(params.template_run_dir).exists()) {
        error "Template directory ${params.template_run_dir} does not exist!"
    }

    existing_params = Channel.fromPath("${params.template_run_dir}/params/*")

    matrices = Channel.fromPath("${params.outdir}/index+matrices.anndata.h5ad")
        | extract_from_anndata 

    out = normalizeMatrix(matrices, existing_params)
}

// De-novo normalization
workflow {
    matrices = Channel.fromPath("${params.outdir}/index+matrices.anndata.h5ad")
        | extract_from_anndata

    out = normalizeMatrix(anndata, Channel.empty())
}
