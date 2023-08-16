#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process split_in_chunks {

    label "highmem"
    conda params.conda

    input:
        path matrix
    
    output:
        path "${prefix}.*"
    
    script:
    prefix = "chunk"
    """
    python3 $moduleDir/bin/split_in_chunks.py \
        ${params.chunksize} \
        ${matrix} \
        ${prefix}
    """
}

process variance_partition {

    conda params.conda
    tag "${chunk_index}"

    input:
        tuple val(chunk_index), path(chunk), path(sample_names)
    
    output:
        path name
    
    script:
    name = "${chunk_index}.variance_partition.tsv"
    """
    Rscript $moduleDir/bin/variance_partition.R \
        ${params.metadata} \
        ${chunk} \
        ${sample_names} \
        var_partition.tsv
    
    cat ${params.filtered_masterlist} \
        | grep -v '#' \
        | sed -n '${chunk_index},${chunk_index + params.chunksize - 1}' \
        | paste - var_partition.tsv > ${name}
    """
}


workflow {
    params.chunksize = 500
    params.samples_order = "$launchDir/${params.outdir}/samples_order.txt"
    params.filtered_masterlist = "$launchDir/${params.outdir}/ masterlist.filtered.bed"
    Channel.fromPath("$launchDir/${params.outdir}/deseq.normalized.sf.vst.npy")
        | split_in_chunks
        | flatten()
        | map(it -> tuple(it.extension.toInteger(), it))
        | variance_partition
        | collectFile(
            name: "materlist.vp_annotated.bed",
            storeDir: params.outdir
        )
}
