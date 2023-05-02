#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.hyper_params_list = "/home/sboytsov/NMF/nmf_hyperparams.tsv"

weights_file_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/sample_weights_annotation_ontology.tsv"
params.matrix_path = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/deseq.normalized.sf.vst.npy"


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

workflow fitNMF {
    take:
        // ID,
        // peaks_id, peaks_params,
        // encoder_id, encoder_params,
        // clustering_alg, clustering_params
        hyperparams 
    main:
        fit_nmf_params = hyperparams
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
            row.n_components,
            row.method))
        | fitNMF
        | visualizeNMF
}

