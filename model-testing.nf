#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process subset_peaks {
    
}

process fit_vae {
	tag "${indiv_id}"
	conda params.conda
    label "gpu"

	input:
		tuple val(vae_params), path(peaks_matrix)

	output:
		tuple val(vae_params), file("${prefix}*")

	script:
    prefix = "vae."
	"""
    #FIXME
	"""
}

process clustering {
    conda params.conda
    tag "${id}"

    input:
        tuple val(clust_alg), val(clust_params), val(id), path(embedding)

    output:
        tuple val(id), path("${prefix}*")
    
    script:
    prefix = "${id}.clustering."
    switch (clust_alg) {
        case "k-means": 
            """

            """
            break;
        case "hierarhical":
            """
            
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

process calc_metrics {
    tag "${id}"
    conda params.conda

    input:
        tuple val(id), path(clustering)
    
    output:
        tuple val(id), path(name)
    script:
    name = "metrics.tsv"
    """
    
    """
}

workflow fitModels {
    take:
        hyperparams
    main:
        peaks = hyperparams 
            | map(it -> it[0])
            | unique
            | subset_peaks
        
        embedding = peaks
            | join(hyperparams)
            | map(it -> tuple(it[1], it[2]))
            | unique
            | fit_vae
        
        out = hyperparams 
            | map(it -> tuple(*it[1..(it.size()-1)]))
            | join(embedding)
            | map(iter -> tuple(*iter[1..(iter.size()-1)]))
            | clustering
            | calc_metrics
    emit:
        out
}


workflow {
    Channel.fromPath("${params.meta_params}")
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.peaks_params, row.encoder_params, row.clust_alg, row.clust_params, row.id))
        | fitModels
}

