#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./normalize_signal"




params.params_list = "/home/sboytsov/NMF/nmf_hyperparams.tsv"
params.samples_order_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/indivs_order.txt"
params.clustering_meta = "/home/sboytsov/poster_clustering/2902_cluster_meta_0303.tsv"



process fit_nmf {
	tag "${prefix}"
	conda "/home/sabramov/miniconda3/envs/jupyterlab"
    publishDir "${params.outdir}/nmf/${prefix}"
    memory { 400.GB * task.attempt }

	input:
		tuple val(prefix), val(n_components), path(matrix_path), val(weights_path), val(peaks_mask), val(samples_mask)

	output:
        tuple val(prefix), val(n_components), path(matrix_path), path("${prefix}.W.npy"), path("${prefix}.H.npy"), path("${prefix}.non_zero_peaks_mask.txt"), path("${prefix}.samples_mask.txt")

	script:
	"""
    python3 $moduleDir/bin/post_processing/perform_NMF.py \
        ${matrix_path} \
        ${prefix} \
        ${n_components} \
        ${non_required_arg(weights_path, '--samples_weights')} \
        ${non_required_arg(samples_mask, '--samples_mask')} \
        ${non_required_arg(peaks_mask, '--peaks_mask')}
	"""
}

process visualize_nmf {
	tag "${prefix}"
	conda "/home/sabramov/miniconda3/envs/jupyterlab"
    publishDir "${params.outdir}/nmf/${prefix}"
    memory { 300.GB * task.attempt }

	input:
        tuple val(prefix), val(n_components), path(binary_matrix), path(W), path(H), path(peaks_mask), path(samples_mask)

	output:
        tuple val(prefix), path("*.pdf")

	script:
	"""
    echo 1
    python3 $moduleDir/bin/post_processing/visualize_nmf.py \
        ${binary_matrix} \
        ${W} \
        ${H} \
        ${params.samples_file} \
        ${n_components} \
        --peaks_mask ${peaks_mask} \
        --samples_mask ${samples_mask} \
        --outpath ./
	"""
}

workflow runNMF {
    take:
        hyperparams 
    main:
        out = hyperparams
            | fit_nmf
            | visualize_nmf
    emit:
        out
}

// nextflow run ~/projects/SuperIndex/nf-index/nmf.nf -profile Altius -resume
workflow {
    Channel.fromPath(params.params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            "${row.prefix}.${row.n_components}",
            row.n_components,
            file(row.matrix_path),
            row?.samples_weights,
            row?.peaks_mask,
            row?.samples_mask
            ))
        | runNMF
}

// Entry for visuzizations only
workflow visualize {
    params.nmf_results_path = "/net/seq/data2/projects/sabramov/SuperIndex/NMF0508/output/nmf_results/"
    Channel.fromPath(params.params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            row.prefix,
            row.n_components,
            row?.samples_mask,
            file("${params.nmf_results_path}/${row.prefix}.${row.n_components}*")
            ))
        | groupTuple(by:[0,1])
        | visualize_nmf
}