#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.params_list = "/home/sboytsov/NMF/nmf_hyperparams.tsv"

params.weights_file_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/sample_weights_annotation_ontology.tsv"
params.matrix_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/deseq.normalized.sf.vst.npy"


process fit_nmf {
	tag "${n_components}:${method}"
	conda "/home/sabramov/miniconda3/envs/jupyter2"
    publishDir "${params.outdir}/nmf_results"
    memory { 400.GB * task.attempt }

	input:
		tuple val(n_components), val(method)

	output:
        tuple val(n_components), val(method), path("${method}*")

	script:
	"""
    python3 $moduleDir/bin/perform_NMF.py \
        ${params.weights_file_path} \
        ${params.matrix_path} \
        ./ \
        ${method} \
        ${n_components}
	"""
}

process visualize_nmf {
	tag "${vae_id}:${peaks_id}"
	conda params.conda
    publishDir "${params.outdir}/visualize"

	input:
		tuple val(n_components), val(method)

	output:
        tuple val(n_components), val(method), path("${prefix}*")

	script:
    prefix = "${method}*"
	"""
    python3 $moduleDir/bin/perfrom_NMF.py \
        ${params.weights_file_path} \
        ${params.matrix_path} \
        ./ \
        ${method} \
        ${n_components}
	"""
}

workflow runNMF {
    take:
        // ID,
        // peaks_id, peaks_params,
        // encoder_id, encoder_params,
        // clustering_alg, clustering_params
        hyperparams 
    main:
        out = fit_nmf(hyperparams) // | visualize_nmf
    emit:
        out
}

workflow visualize {
    Channel.fromPath('/net/seq/data2/projects/sabramov/SuperIndex/NMF0501/output/nmf_results/*')
        | map(it -> tuple(it.name.split('.', 0), it.name.split('.', 0), it))
        | view()
        | groupTuple(by:[0,1])
        | visualize_nmf
}


workflow {
    Channel.fromPath(params.params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.n_components, row.method))
        | runNMF
}

