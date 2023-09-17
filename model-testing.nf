#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.hyper_params_list = "/home/sabramov/projects/SuperIndex/hyperparams_clustering+ids.tsv"
params.gpu_conda = "/home/sabramov/miniconda3/envs/tensorflow"

params.meta = "/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_22-11-28/metadata/ENCODE4plus_master_metadata_filtered.tsv"

params.normalized_matrices = "$launchDir/${params.outdir}/normalized_matrices.hdf5"


process subset_peaks {
    conda params.conda
    tag "${id}"
    publishDir "${params.outdir}/peaks", pattern: "${name}"
    publishDir "${params.outdir}/peaks_masks", pattern: "${prefix}.mask.txt"
    memory 350.GB
    input:
		tuple val(id), val(peaks_params)

	output:
		tuple val(id), path(name), emit: matrix
        path "${prefix}.mask.txt", emit: mask

    
    script:
    prefix = "${id}.peaks"
    name = "${prefix}.npy"
    """
    echo -e '${peaks_params}' > params.json
    python3 $moduleDir/bin/subset_peaks.py \
        params.json \
        ${params.normalized_matrix} \
        ${prefix}
    """
}

process fit_vae {
	tag "${vae_id}:${peaks_id}"
	conda params.gpu_conda
    label "gpu"
    publishDir "${params.outdir}/vae", pattern: "${name}"
    publishDir "${params.outdir}/vae_models", pattern: "${prefix}_*"
    scratch true

	input:
		tuple val(peaks_id), val(vae_id), val(vae_params), path(peaks_matrix)

	output:
        tuple val(peaks_id), val(vae_id), path(name), emit: emb
		path "${prefix}_*", emit: all_data

	script:
    prefix = "${vae_id}.${peaks_id}"
    name = "${prefix}.embedding.npy"
	"""
    echo '${vae_params}' > params.json
    python3 $moduleDir/bin/fit_auto_encoder.py \
        params.json \
        ${peaks_matrix} \
        ${prefix}
	"""
}

process clustering {
    conda params.conda
    tag "${id}"
    publishDir "${params.outdir}/clustering", pattern: "${name}"
    publishDir "${params.outdir}/clustering_data", pattern: "${id}.[0-9]*"


    input:
        tuple val(peaks_id), val(encoder_id), val(id), val(clust_alg), val(clust_params), path(embedding)

    output:
        tuple val(id), path("${name}"), emit: metrics
        tuple val(id), path("${id}*"), emit: all_data
    
    script:
    name = "${id}.clustering.metrics.tsv"
    switch (clust_alg) {
        case "k_means": 
            """
            echo '${clust_params}' > params.json
            python3 $moduleDir/bin/k-means.py params.json ${embedding} ${id}
            """
            break;
        case "hierarchical":
            """
            echo '${clust_params}' > params.json
            python3 $moduleDir/bin/aggloclustering.py \
                params.json \
                ${embedding} \
                ${params.meta} \
                ${params.normalized_matrices} \
                ${id}
            """
            break;
        case "community":
            """
            exit 1
            """
            break;
        default: 
            error "Clustering with ${clust_alg} is not implemented."
            break;
    }

}

workflow fitModels {
    take:
        // ID,
        // peaks_id, peaks_params,
        // encoder_id, encoder_params,
        // clustering_alg, clustering_params
        hyperparams 
    main:
        subset_peaks_params = hyperparams
            | map(it -> tuple(it[1], it[2]))
            | unique()
        peaks = subset_peaks(subset_peaks_params).matrix // peaks_id, peaks_subset
        
        embedding = hyperparams
            | map(it -> tuple(it[1], it[3], it[4])) // peaks_id, encoder_id, encoder_params
            | unique()
            | combine(peaks, by: 0) // peaks_id, encoder_id, encoder_params, peaks_subset
            | fit_vae // peaks_id, encoder_id, embedding

        out = hyperparams 
            | map(it -> tuple(it[1], it[3], it[0], it[5], it[6])) //  peaks_id, encoder_id, clustering_alg, clustering_params
            | combine(embedding.emb, by: [0, 1]) //  peaks_id, encoder_id, ID, clustering_alg, clustering_params, embedding
            | clustering
    emit:
        out.metrics
}


workflow {
    Channel.fromPath(params.hyper_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.id,
            row.peak_id,
            row.peaks_params,
            row.encoder_id,
            row.encoder_params,
            row.clust_alg,
            row.clust_params))
        | fitModels
        | map(it -> it[1])
        | collectFile(name: 'all.metrics.tsv', 
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)
}


process select_peaks {

    publishDir "${params.outdir}/${prefix}"
    conda params.conda
    tag "${prefix}"

    input:
        tuple val(prefix), val(peaks_params), path(signal_matrix), path(binary_matrix), path(peaks_meta), path(samples_meta)
    
    output:
        path("${prefix}*")

    script:
    """
    echo -e '${peaks_params}' > params.json
    python3 $moduleDir/bin/select_peaks.py \
        params.json \
        ${peaks_meta} \
        ${samples_meta} \
        ${signal_matrix} \
        ${binary_matrix} \
        ${prefix}
    """

}

workflow selectPeaks {
    Channel.fromPath(params.hyper_params_list)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
            row.prefix,
            row.params,
            row.signal_matrix,
            row.binary_matrix,
            row.peaks_meta,
            row.samples_meta))
        | select_peaks
}

