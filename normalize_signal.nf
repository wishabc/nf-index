#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { variancePartition } from "./variance_partition"
include { add_normalized_matrices_to_anndata } from "./convert_to_anndata"

params.conda = "$moduleDir/environment.yml"
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

		normalization = normalize_matrix(
            matrices.map(it -> tuple(it[0], it[1])),
            lowess_params
        )

		out = deseq2(normalization.scale_factors, deseq_params).matrix
            | combine(deseq2.out.model_params)
            | combine(normalization.log_diffs)
            | combine(sf)
            | combine(normalization.model_params)

	emit:
		out // deseq2_matrix, deseq2_model_params, log_diffs, scale_factors, lowess_params1, lowess_params2
}

// Re-use existing model from other run (e.g. different samples)
workflow existingModel {
    params.template_run_dir = "/path/to/previous/run"
    if (!file(params.template_run_dir).exists()) {
        error "Template directory ${params.template_run_dir} does not exist!"
    }

    existing_params = Channel.fromPath("${params.template_run_dir}/params/*")

    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)

    out = normalizeMatrix(matrices, existing_params)
        | map(it -> tuple([it[0], it[2], it[3]], [it[1], it[4], it[5]]))
 
    add_normalized_matrices_to_anndata(anndata, out)
}

// De-novo normalization
workflow {
    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)

    out = normalizeMatrix(matrices, Channel.empty())

    add_normalized_matrices_to_anndata(anndata, out)
}
