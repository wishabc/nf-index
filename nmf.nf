#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./normalize_signal"


process fit_nmf {
	tag "${prefix}"
	conda params.conda
    publishDir "${params.outdir}/nmf/${prefix}", pattern: "${prefix}.*"
    label "highmem"

	input:
		tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), val(samples_weights), val(peaks_mask), val(samples_mask), val(peaks_weights), val(extra_params)

	output:
        tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), path("${prefix}.W.npy"), path("${prefix}.H.npy"), path("${prefix}.non_zero_peaks_mask.txt"), path("${prefix}.samples_mask.txt")

	script:
	"""
    python3 $moduleDir/bin/nmf/perform_NMF.py \
        ${n_components} \
        ${prefix} \
        --matrix ${matrix_path} \
        --sample_names ${sample_names} \
        --dhs_meta ${dhs_meta} \
        ${non_required_arg(samples_weights, '--samples_weights')} \
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
        tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), path(W), path(H), path(peaks_mask), path(samples_mask)

	output:
        tuple val(prefix), path("*.pdf")

	script:
	"""
    python3 $moduleDir/bin/nmf/visualize_nmf.py \
        ${n_components} \
        ${W} \
        ${H} \
        --outpath ./ \
        --matrix ${matrix_path} \
        --sample_names ${sample_names} \
        --dhs_meta ${dhs_meta} \
        --samples_metadata ${params.samples_file} \
        --peaks_mask ${peaks_mask} \
        --samples_mask ${samples_mask} \
        ${non_required_arg(params.dhs_annotations, '--dhs_annotations')}
	"""
}

process add_metadata {
    conda params.conda
    publishDir "${params.outdir}"

    output:
        path name

    script:
    name = "${file(params.nmf_params_list).baseName}+matrices.tsv"
    """
    python3 $moduleDir/bin/nmf/add_metadata.py \
        ${params.nmf_params_list} \
        ${params.outdir}/nmf \
        ${name}
    """
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
