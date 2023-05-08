#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.params_list = "/home/sboytsov/NMF/nmf_hyperparams.tsv"


params.samples_order_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/indivs_order.txt"
params.clustering_meta = "/home/sboytsov/poster_clustering/2902_cluster_meta_0303.tsv"

def non_required_arg(value, key) {
    return value ? "${key} ${value}": ""
}

process fit_nmf {
	tag "${prefix}:${n_components}"
	conda "/home/sabramov/miniconda3/envs/jupyter2"
    publishDir "${params.outdir}/nmf_results"
    memory { 400.GB * task.attempt }

	input:
		tuple val(n_components), val(fname), path(matrix_path), val(weights_path), val(peaks_mask), val(samples_mask)

	output:
        tuple val(prefix), path("${prefix}*")

	script:
    weights = weights_path ? "--sampels_weights ${weights_path}": ""
    prefix = "${fname}.${n_components}"
	"""
    python3 $moduleDir/bin/perform_NMF.py \
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
	conda params.conda
    errorStrategy 'ignore'
    publishDir "${params.outdir}/figures"

	input:
		tuple val(prefix), path(nmf_results)

	output:
        tuple val(prefix), path("*.pdf")

	script:
	"""
    python3 $moduleDir/bin/visualize_nmf.py \
        ${params.clustering_meta} \
        ${params.samples_order_path} \
        ${prefix} \
        ${n_components}
	"""
}

workflow runNMF {
    take:
        hyperparams 
    main:
        out = fit_nmf(hyperparams) // | visualize_nmf
    emit:
        out
}

workflow visualize {
    data = Channel.fromPath('/net/seq/data2/projects/sabramov/SuperIndex/NMF0508/output/nmf_results/*')
        | map(it -> tuple(it.name.split('\\.')[2], it.name.split('\\.')[0], it))
        | groupTuple(by:[0,1])
        | visualize_nmf
}


workflow {
    Channel.fromPath(params.params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            row.n_components,
            row.prefix,
            file(row.matrix_path),
            row?.weights_path,
            row?.peaks_mask,
            row?.samples_mask
            ))
        | runNMF
}

