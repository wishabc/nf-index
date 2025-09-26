#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { variancePartition } from "./variance_partition"
include { add_normalized_matrices_to_anndata } from "./converters"

params.conda = "$moduleDir/environment.yml"
params.sample_weights = ""


def non_required_arg(value, key) {
    return value ? "${key} ${value}": ""
}


process extract_from_anndata {
    conda params.conda
    label "bigmem"

    input:
        val anndata
    
    output:
        tuple path("binary.matrix.npy"), path("counts.matrix.npy"), path(samples_order), path(masterlist), path(samples_meta)

    script:
    samples_order = "samples_order.txt"
    masterlist = "masterlist.no_header.bed"
    samples_meta = "samples_meta.txt"
    """
    python3 $moduleDir/bin/convert_to_anndata/extract_from_anndata.py \
        ${anndata} \
        ${masterlist} \
        ${samples_order} \
        ${samples_meta} \
        --extra_layers binary counts \
        --dhs_mask_name autosomal_dhs
    """
}


process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}", pattern: "${prefix}*.npy"
    publishDir "${params.outdir}/params", pattern: "${prefix}.lowess_params*"
    publishDir "${params.outdir}/qc", pattern: "${prefix}*.pdf"

	input:
        tuple path(peaks_matrix), path(signal_matrix)
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
	n = norm_params.name != "params/empty.params" ? file(norm_params[0]).baseName : ""
	normalization_params = n ? "--model_params params/${n}" : ""
	"""
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1

	normalize-matrix ${signal_matrix} \
        ${save_dir} \
        --prefix ${pref} \
		--jobs ${task.cpus} \
        --normalization_type lowess \
		${normalization_params}

    qc-normalization ${save_dir} \
        --prefix ${pref} \
        --n_samples 10 \
        --peak_matrix ${peaks_matrix}
	"""
}

process prepare_data_for_deseq2 {
	conda params.conda
	publishDir "${params.outdir}/params", pattern: "${prefix}.params.RDS"
    publishDir "${params.outdir}/deseq2", pattern: "${prefix}.dds.RDS"
	label "bigmem"

	input:
		tuple path(scale_factors), path(signal_matrix), path(samples_order)
		path norm_params, stageAs: "params/*"

	output:
		tuple val(prefix), path("${prefix}.dds.RDS"), emit: data
        path "${prefix}.params.RDS", emit: model_params

	script:
	prefix = "deseq_normalized.only_autosomes.filtered.sf.vst"
	normalization_params = norm_params.name != "params/empty.params" ? norm_params : ""
	"""
	Rscript $moduleDir/bin/prepare_data_for_deseq.R \
		${signal_matrix} \
		${scale_factors} \
		${samples_order} \
		${prefix} \
		${normalization_params}
	"""
}

process deseq2_vst {
	conda params.conda
	publishDir "${params.outdir}/normalization"
	label "bigmem"

	input:
		tuple val(prefix), path(dataset)

	output:
		path name

	script:
    name = "${prefix}.npy"
	"""
	Rscript $moduleDir/bin/deseq2_vst.R \
		${dataset} \
		${name} 
	"""
}


workflow normalizeMatrix {
	take:
		matrices  // binary_matrix, count_matrix, samples_order, masterlist, samples_file
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

        dat = normalization.scale_factors
            | combine(matrices)
            | map(it -> tuple(it[0], it[2], it[3]))
        
        normalized_matrix = prepare_data_for_deseq2(dat, deseq_params).data
            | deseq2_vst

        vp = normalized_matrix
            | combine(matrices)
            | map(it -> tuple(it[0], it[4], it[5], params.formula))
            | variancePartition

		out = normalized_matrix
            | combine(prepare_data_for_deseq2.out.model_params)
            | combine(normalization.model_params)
            | combine(vp)
            | map(it -> tuple([it[0]], [it[1], it[2], it[3]], it[4], it[5]))

	emit:
		out // deseq2_matrix, [model_params], vp_annotated_masterlist, formula
}

process extract_normalization_params_from_template {
    conda params.conda
    label "medmem"

    input:
        val anndata

    output:
        path("params/*")

    script:
    """
    mkdir params
    python3 $moduleDir/bin/convert_to_anndata/extract_normalization_params.py \
        ${anndata} \
        params
    """
}

// Re-use existing model from other run (e.g. different samples)
workflow existingModel {
    params.template_anndata = "/path/to/previous/anndata_with_params"
    if (!file(params.template_anndata).exists()) {
        error "Template anndata: ${params.template_anndata} does not exist!"
    } else {
        print "Using existing model from ${params.template_anndata}"
    }

    existing_params = Channel.of(params.template_anndata)
        | extract_normalization_params_from_template
        | flatten()

    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)

    out = normalizeMatrix(matrices, existing_params)
        
    add_normalized_matrices_to_anndata(anndata, out)
}

// De-novo normalization
workflow {
    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)

    out = normalizeMatrix(matrices, Channel.empty())

    add_normalized_matrices_to_anndata(anndata, out)
}


process differential_deseq {

    memory { 30.GB * task.attempt * task.attempt }
    conda "/home/sabramov/miniconda3/envs/r-jupyter/"
    label "medmem"

    input:
        val chunk

    output:
        path "*${suffix}.tsv"

    script:
    suffix = "deseq_res.${chunk}"
    """
    Rscript $moduleDir/bin/differential_deseq.R \
        ${suffix} \
        ${params.dds} \
        ${params.total_chunks} \
        ${chunk}
    """
}


workflow diffDeseq {
    params.dds = "/net/seq/data2/projects/sabramov/SuperIndex/hotspot3/w_babachi_new.v23/index/filled_1pr/output/normalization/lowess/ready_for_deseq.reciprocal.mouse+human.normalized.only_autosomes.filtered.no_q.dds.RDS"
    params.total_chunks = 500
    Channel.of(1..params.total_chunks)
        | differential_deseq
        | flatten()
        | collectFile(
            storeDir: "${params.outdir}",
            skip: 1,
            keepHeader: true
        ) {
            it -> [ "${it.simpleName}.tsv", it ]
        }
}
