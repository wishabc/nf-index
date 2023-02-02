#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TODO: move to main.nf
process filter_singletons {
    conda params.conda
    publishDir "${params.outdir}/index", pattern: "${name}"
    scratch true

    output:
        path name

    script:
    name = "singletons_mask.txt"
    """
    cat ${params.index_file} | grep -v chrX | grep -v chrY > filtered_index.bed
    bedmap --indicator --sweep-all --bp-ovr 1 filtered_index.bed \
        ${params.encode_blacklist_regions} > blacklisted_mask.txt
    
    python3 $moduleDir/bin/filter_index.py \
        filtered_index.bed \
        blacklisted_mask.txt \
        ${name}
    """
}


process subset_peaks {
    conda params.conda
    tag "${id}"
    publishDir "${params.outdir}/matrix", pattern: "${name}"
    memory 350.GB
    input:
		tuple val(id), val(peaks_params)
        path singletons_mask

	output:
		tuple val(peaks_params), path(name)
    
    script:
    name = "${id}.peaks.npy"
    """
    echo -e '${peaks_params}' > params.json
    python3 $moduleDir/bin/subset_peaks.py \
        params.json \
        ${params.normalized_matrix} \
        ${singletons_mask} \
        ${name}
    """
}

process fit_vae {
	tag "${id}"
	conda params.gpu_conda
    label "gpu"
    publishDir "${params.outdir}/vae", pattern: "${name}"
    publishDir "${params.outdir}/vae_models", pattern: "${id}_*"
    scratch true

	input:
		tuple val(id), val(vae_params), path(peaks_matrix)

	output:
        tuple val(vae_params), path(name), emit: emb
		path "${id}_*", emit: all_data

	script:
    name = "${id}.embedding.npy"
	"""
    echo '${vae_params}' > params.json
    python3 $moduleDir/bin/fit_auto_encoder.py \
        params.json \
        ${peaks_matrix} \
        ${id}
	"""
}

process clustering {
    conda params.conda
    tag "${id}"
    publishDir "${params.outdir}/clustering", pattern: "${prefix}.metrics.tsv"
    publishDir "${params.outdir}/clustering_data", pattern: "${prefix}.[0-9]*"

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
            echo '${clust_params}' > params.json
            python3 $moduleDir/bin/hierarchical.py params.json ${embedding} ${prefix}
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
        hyperparams // ID, peaks_params, encoder_params, encoder_params, clustering_alg, clustering_params
    main:
        out_mask = filter_singletons()
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
    params.normalized_matrix = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/deseq.normalized.vst.txt.npy"
    params.meta_params = "/home/sabramov/projects/SuperIndex/hyperparams_clustering.tsv"
    Channel.fromPath(params.meta_params)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.id, row.peaks_params,
            row.encoder_params, row.clust_alg, row.clust_params))
        | fitModels
}

