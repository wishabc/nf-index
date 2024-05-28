#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./normalize_signal"


process fit_nmf {
	tag "${prefix}"
	conda "/home/sabramov/miniconda3/envs/jupyterlab"
    publishDir "${params.outdir}/nmf/${prefix}", pattern: "${prefix}.*"
    label "highmem"

	input:
		tuple val(prefix), val(n_components), path(matrix_path), path(sample_names), path(dhs_meta), val(weights_path), val(peaks_mask), val(samples_mask), val(peaks_weights)

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
        ${non_required_arg(peaks_weights, '--peaks_weights')}
	"""
}

process visualize_nmf {
	tag "${prefix}"
	conda "/home/sabramov/miniconda3/envs/jupyterlab"
    publishDir "${params.outdir}/nmf/${prefix}"
    label "highmem"

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

process top_samples_track {

    scratch true
    conda params.conda
    tag "${component}"
    publishDir "${params.outdir}/top_samples"

    input:
        tuple val(component), path(density_bw, stageAs: "?/*")
    
    output:
        tuple path(name), path("${component}.top_samples.bg")
    
    script:
    name = "${component}.top_samples.bw"
    bg = "${component}.top_samples.bg"
    """
    wiggletools write_bg ${bg} mean ${density_bw}
    bedGraphToBigWig "${bg}" "${params.chrom_sizes}" "${name}"
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
            row?.peaks_weights
            ))
        | runNMF
}

// Entry for visuzizations only
workflow visualize {
    params.nmf_results_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-index0415/matrices/downsampled_no_cancer/output/nmf/"
    Channel.fromPath(params.nmf_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            "${row.prefix}.${row.n_components}",
            row.n_components,
            file(row.matrix_path),
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
        | visualize_nmf
}


workflow topSamples {
    meta = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.normalized_density_bw)))
    Channel.fromPath(params.top_components_ids)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, row.component))
        | join(meta)
        | map(it -> tuple(it[1], it[2]))
        | groupTuple()
        | top_samples_track
}
