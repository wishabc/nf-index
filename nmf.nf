#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./normalize_signal"


process fit_nmf {
	tag "${prefix}"
	conda params.conda
    publishDir "${params.outdir}/nmf/${prefix}", pattern: "${prefix}.*"
    label "highmem"

	input:
		tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), val(weights_path), val(peaks_mask), val(samples_mask), val(peaks_weights), val(extra_params)

	output:
        tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), path("${prefix}.W.npy"), path("${prefix}.H.npy"), path("${prefix}.non_zero_peaks_mask.txt"), path("${prefix}.samples_mask.txt")

	script:
	"""
    python3 $moduleDir/bin/post_processing/perform_NMF.py \
        ${matrix_path} \
        ${sample_names} \
        ${dhs_meta} \
        ${prefix} \
        ${n_components} \
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
        tuple val(prefix), val(n_components), path(binary_matrix), path(sample_names), path(masterlist), path(W), path(H), path(peaks_mask), path(samples_mask)

	output:
        tuple val(prefix), path("*.pdf")

	script:
	"""
    python3 $moduleDir/bin/post_processing/visualize_nmf.py \
        ${binary_matrix} \
        ${sample_names} \
        ${W} \
        ${H} \
        ${params.samples_file} \
        ${masterlist} \
        ${n_components} \
        --peaks_mask ${peaks_mask} \
        --samples_mask ${samples_mask} \
        --outpath ./ \
        ${non_required_arg(params.dhs_annotations, '--dhs_annotations')}
	"""
}

process add_metadata {
    conda params.conda
    publishDir "${params.outdir}"

    output:
        tuple path(name)

    script:
    name = "nmf_meta.tsv"
    """
    python3 $moduleDir/bin/post_processing/add_metadata.py \
        ${params.samples_file} \
        ${params.outdir}/nmf \
        ${name}
    """
}


workflow runNMF {
    take:
        hyperparams 
    main:
        out = hyperparams
            | fit_nmf
            | visualize_nmf
        add_metadata()
    emit:
        out
}

// nextflow run ~/projects/SuperIndex/nf-index/nmf.nf -profile Altius -resume
workflow {
    Channel.fromPath(params.nmf_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            "${row.prefix}.${row.n_components}",
            row.n_components,
            file(row.matrix_path),
            file(row.sample_names),
            file(row.dhs_meta),
            row?.samples_weights,
            row?.peaks_mask,
            row?.samples_mask,
            row?.peaks_weights,
            row?.extra_params
            ))
        | runNMF
}

// Entry for visuzizations only
workflow visualize {
    Channel.fromPath(params.nmf_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            "${row.prefix}.${row.n_components}",
            row.n_components,
            file(row.matrix_path),
            file(row.sample_names),
            file(row.dhs_meta)
            ))
        | map( 
            it -> tuple(
                *it[0..(it.size()-1)],
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.W.npy"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.H.npy"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.non_zero_peaks_mask.txt"),
                file("${params.nmf_results_path}/${it[0]}/${it[0]}.samples_mask.txt"),
            )
        )
        | filter { it[5].exists() }
        | visualize_nmf
}
