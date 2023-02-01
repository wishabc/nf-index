#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TODO: move to main.nf
process filter_singletons {
    conda params.conda
    publishDir "${params.outdir}/index", pattern: name
    scratch true

    output:
        path name

    script:
    name = "singletons_mask.txt"
    """
    bedmap --indicator ${params.index_file} ${params.encode_blacklist_regions} > blacklisted_mask.txt
    python3 $moduleDir/bin/filter_index.py ${params.index_file} blacklisted_mask.txt ${name}
    """
}


process subset_peaks {
    conda params.conda
    tag "${id}"
    publishDir "${params.outdir}/matrix", pattern: name
    
    input:
		tuple val(id), val(peaks_params)
        path singletons_mask

	output:
		tuple val(peaks_params), path(name)
    
    script:
    name = "${id}.peaks.npy"
    """
    echo "${peaks_params}" > params.json
    python3 $moduleDir/bin/subset_peaks.py params.json ${params.normalized_matrix} ${singletons_mask} ${name}
    """
}

process fit_vae {
	tag "${id}"
	conda params.conda
    label "gpu"
    publishDir "${params.outdir}/vae", pattern: "${id}.npy"
    publishDir "${params.outdir}/vae_models", pattern: "${id}.model.*"

	input:
		tuple val(id), val(vae_params), path(peaks_matrix)

	output:
        tuple val(vae_params), path("${id}.model.npy"), emit: emb
		path "${prefix}*", emit: all_data

	script:
	"""
    # TODO save model
    echo "${vae_params}" > params.json
    python3 $moduleDir/bin/fit_auto_encoder.py \
        params.json \
        ${peaks_matrix} \
        ${id}
	"""
}

process clustering {
    conda params.conda
    tag "${id}"

    input:
        tuple val(id), val(clust_alg), val(clust_params), path(embedding)

    output:
        tuple val(id), path("${prefix}*")
    
    script:
    prefix = "${id}.clustering."
    switch (clust_alg) {
        case "k-means": 
            """
            echo "${clust_params}" > params.json
            python3 $moduleDir/bin/k-means.py params.json ${embedding} ${prefix}
            """
            break;
        case "hierarchical":
            """
            echo "${clust_params}" > params.json
            python3 $moduleDir/bin/hierarchical.py params.json ${embedding} ${prefix}
            """
            break;
        case "community":
            """
            
            """
            break;
        default: 
            error "Clustering with ${clust_alg} is not implemented."
            break;
    }

}

workflow fitModels {
    take:
        hyperparams // ID, peaks_params, encoder_params, encoder_params, clustering_alg, clustering_params
    main:
        out_mask = filter_singletons()
        params.normalized_matrix = ""
        subset_peaks_params = hyperparams 
            | map(it -> tuple(it[0], it[1]))
            | unique { it[1] }
        peaks = subset_peaks(subset_peaks_params, out_mask) // peaks_params, peaks_subset
        
        embedding = hyperparams
            | map(it -> tuple(it[1], it[0], it[2]))
            | join(peaks) // peaks_params, ID, encoder_params, peaks_subset
            | map(it -> tuple(*it[1..(it.size()-1)])) // ID, encoder_params, peaks_subset
            | unique { it[1] }
            | fit_vae // encoder_params, embedding

        out = hyperparams 
            | map(it -> tuple(it[3], it[0], it[4], it[5])) //  encoder_params, ID,  clustering_alg, clustering_params
            | join(embedding.emb) //  encoder_params, ID,  clustering_alg, clustering_params, embedding
            | map(it -> tuple(*it[1..(it.size()-1)])) // ID,  clustering_alg, clustering_params, embedding
            | clustering
    emit:
        out
}


workflow {
    Channel.fromPath("${params.meta_params}")
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.id, row.peaks_params,
            row.encoder_params, row.clust_alg, row.clust_params))
        | fitModels
}

