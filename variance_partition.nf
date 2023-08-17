#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process variance_partition {

    conda params.conda
    tag "${start_index}"
    publishDir "${params.outdir}/variance_partition"

    input:
        tuple val(start_index), path(masterlist), path(h5file)
    
    output:
        path name
    
    script:
    name = "${start_index}.variance_partition.tsv"
    end_index = start_index + params.chunk_size - 1
    """
    Rscript $moduleDir/bin/variance_partition.R \
        ${params.samples_file} \
        ${start_index} \
        ${params.chunk_size} \
        ${h5file} \
        ${masterlist} \
        ${name}
    """
}

workflow variancePartition {
    take:
        data // masterlist, h5file
    main:
        total_dhs = data.first().map(it -> it[0]).countLines()
        out = Channel.of(1..total_dhs)
            | collate(params.chunk_size)
            | map(it -> it[0])
            | combine(data)
            | variance_partition
            | collectFile(
                name: "masterlist.vp_annotated.bed",
                storeDir: params.outdir,
                keepHeader: true,
                sort: true,
                skip: 1
            )
    emit:
        out  
}


workflow {
    params.chunk_size = 10000
    params.h5file = "$launchDir/${params.outdir}/matrices.h5"
    params.filtered_masterlist = "$launchDir/${params.outdir}/masterlist.filtered.bed"
    Channel.fromPath(params.filtered_masterlist)
        | combine(
            Channel.fromPath(params.h5file)
        )
        | variancePartition
}
