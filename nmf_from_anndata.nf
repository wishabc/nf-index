#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./normalize_signal"
include { add_metadata } from "./nmf"


process fit_nmf {
	tag "${prefix}"
	conda params.conda
    publishDir "${params.outdir}/nmf/${prefix}", pattern: "${prefix}.*"
    label "highmem"

	input:
		tuple val(prefix), val(n_components), path(anndata_path), val(weights_path), val(peaks_mask), val(samples_mask), val(peaks_weights), val(extra_params)

	output:
        tuple val(prefix), val(n_components), path(anndata_path), path("${prefix}.W.npy"), path("${prefix}.H.npy"), path("${prefix}.non_zero_peaks_mask.txt"), path("${prefix}.samples_mask.txt")

	script:
	"""
    python3 $moduleDir/bin/nmf/perform_NMF.py \
        ${n_components} \
        ${prefix} \
        --from_anndata ${anndata_path} \
        --samples_mask_column reference_sample \
        --dhs_mask_column autosomal_dhs \
        ${non_required_arg(weights_path, '--samples_weights')} \
        ${non_required_arg(samples_mask, '--samples_mask')} \
        ${non_required_arg(peaks_mask, '--peaks_mask')} \
        ${non_required_arg(peaks_weights, '--peaks_weights')} \
        ${non_required_arg(extra_params, '--extra_params')}
	"""
}

process visualize_nmf {
	tag "${prefix}"
	conda params.conda
    publishDir "${params.outdir}/nmf/${prefix}"
    label "highmem"
    errorStrategy 'ignore'

	input:
        tuple val(prefix), val(n_components), path(anndata_path), path(W), path(H), path(peaks_mask), path(samples_mask)

	output:
        tuple val(prefix), path("*.pdf")

	script:
	"""
    python3 $moduleDir/bin/nmf/visualize_nmf.py \
        ${n_components} \
        ${prefix} \
        ${W} \
        ${H} \
        --outpath ./ \
        --from_anndata ${anndata_path} \
        --peaks_mask ${peaks_mask} \
        --samples_mask ${samples_mask} \
        --samples_mask_column reference_sample \
        --dhs_mask_column autosomal_dhs \
	"""
}


// nextflow run ~/projects/SuperIndex/nf-index/nmf.nf -profile Altius -resume
workflow {
    Channel.fromPath(params.nmf_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            row.prefix,
            row.n_components,
            file(row.anndata_path),
            row?.samples_weights,
            row?.peaks_mask,
            row?.samples_mask,
            row?.peaks_weights,
            row?.extra_params
            ))
        | distinct { it[0] }
        | fit_nmf
        | visualize_nmf
    add_metadata()
}


// Entry for visuzizations only
workflow visualize {
    println "Visualizing NMF results from params.nmf_params_list = ${params.nmf_params_list}"
    Channel.fromPath(params.nmf_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            row.prefix,
            row.n_components,
            file(row.anndata_path),
            ))
        | distinct { it[0] }
        | map( 
            it -> tuple(
                *it[0..(it.size()-1)],
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.W.npy"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.H.npy"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.non_zero_peaks_mask.txt"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.samples_mask.txt"),
            )
        )
        | filter { it[3].exists() }
        | visualize_nmf
}
