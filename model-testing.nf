#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


params.meta = "/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_22-11-28/metadata/ENCODE4plus_master_metadata_filtered.tsv"
params.normalized_matrix = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/deseq.normalized.vst.txt.npy"
params.meta_params = "/home/sabramov/projects/SuperIndex/hyperparams_clustering+ids.tsv"

params.indivs_order = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/indivs_order.txt"


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
		tuple val(id), path(name)
    
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
		tuple val(peaks_id), val(id), val(vae_params), path(peaks_matrix)

	output:
        tuple val(peaks_id), val(id), path(name), emit: emb
		path "${id}_*", emit: all_data

	script:
    name = "${id}.embedding.npy"
	"""
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/
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
        tuple val(peaks_id), val(encoder_id), val(id), val(clust_alg), val(clust_params), path(embedding)

    output:
        tuple val(id), path("${prefix}.metrics.tsv"), emit: metrics
        tuple val(id), path("${prefix}*"), emit: all_data
    
    script:
    prefix = "${id}.clustering"
    switch (clust_alg) {
        case "k_means": 
            """
            echo "${clust_params}" > params.json
            python3 $moduleDir/bin/k-means.py params.json ${embedding} ${prefix}
            """
            break;
        case "hierarchical":
            """
            echo '${clust_params}' > params.json
            python3 $moduleDir/bin/aggloclustering.py \
                params.json \
                ${embedding} \
                ${params.meta} \
                ${params.indivs_order} \
                ${prefix}
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
        out_mask = filter_singletons()

        subset_peaks_params = hyperparams
            | map(it -> tuple(it[1], it[2]))
            | unique()
        peaks = subset_peaks(subset_peaks_params, out_mask) // peaks_id, peaks_subset
        
        embedding = hyperparams
            | map(it -> tuple(it[1], it[3], it[4])) // peaks_id, encoder_id, encoder_params
            | unique()
            | combine(peaks, by: 0) // peaks_id, encoder_id, encoder_params, peaks_subset
            | fit_vae // peaks_id, encoder_id, embedding

        out = hyperparams 
            | map(it -> tuple(it[1], it[3], it[0], it[5], it[6], it[7])) //  peaks_id, encoder_id, clustering_alg, clust_id, clustering_params
            | combine(embedding.emb, by: [0, 1]) //  peaks_id, encoder_id, ID, clustering_alg, clustering_params, embedding
            | clustering
    emit:
        out.metrics
}


workflow {
    Channel.fromPath(params.meta_params)
        | splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.id,
            row.peak_id,
            row.peaks_params,
            row.encoder_id,
            row.encoder_params,
            row.clust_alg,
            row.clust_params))
        | fitModels
        | collectFile(name: "all.metrics.tsv", 
            storeDir: "${params.outdir}/all_data",
            skip: 1,
            sort: true,
            keepHeader: true)
}

