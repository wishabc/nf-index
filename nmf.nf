#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.params_list = "/home/sboytsov/NMF/nmf_hyperparams.tsv"

params.weights_file_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/sample_weights_annotation_ontology.tsv"
params.matrix_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/deseq.normalized.sf.vst.npy"
params.sample_order_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/indivs_order.txt"
params.meta_path = "/home/sabramov/projects/SuperIndex/index_clustering_2023-02-08/ENCODE4_altius_index_clustering_metadata_2023-02-08.tsv"
params.cluster_meta_path = "/home/sboytsov/poster_clustering/2902_cluster_meta_0303.tsv"
params.gen_meta_path = "/home/sabramov/projects/ENCODE4/release_0103/genotyping_meta_230206+ids.tsv"

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
    publishDir "${params.outdir}"

	input:
		tuple val(n_components), val(method), path(nmf_results)

	output:
        tuple val(n_components), val(method), path("./figures/*")

	script:
	"""
    python3 $moduleDir/bin/visualize_nmf.py \
        ${params.sample_order_path} \
        ${params.meta_path} \
        ${params.cluster_meta_path} \
        ${params.gen_meta_path} \
        ./ \
        ./figures/ \
        ${method} \
        ${n_components}
	"""
}

workflow runNMF {
    take:
        hyperparams 
    main:
        out = fit_nmf(hyperparams) | visualize_nmf
    emit:
        out
}

workflow visualize {
    Channel.fromPath('/net/seq/data2/projects/sabramov/SuperIndex/NMF0501/output/nmf_results/*')
        | map(it -> tuple(it.name.split('\\.')[0], it.name.split('\\.')[2], it))
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
