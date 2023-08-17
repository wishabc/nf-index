#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process variance_partition {

    conda params.conda
    tag "${start_index}"

    input:
        val start_index
    
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
        ${params.h5file} \
        var_partition.tsv

    tail -n +2 ${params.filtered_masterlist} \
        | sed -n '${start_index},${end_index} p' > masterlist_slice.bed
    
    head -1 ${params.filtered_masterlist} > header_masterlist.bed 
    head -1 var_partition.tsv | paste header_masterlist.bed  - > ${name}
    tail -n +2 var_partition.tsv \
        | paste masterlist_slice.bed - >> ${name}
    
    
    """
}


workflow {
    params.chunk_size = 10000
    params.h5file = "$launchDir/${params.outdir}/matrices.h5"
    
    params.filtered_masterlist = "$launchDir/${params.outdir}/masterlist.filtered.bed"
    total_dhs = file(params.filtered_masterlist).countLines()
    Channel.of(1..total_dhs)
        | collate(params.chunk_size)
        | map(it -> it[0])
        | variance_partition
        | collectFile(
            name: "materlist.vp_annotated.bed",
            storeDir: params.outdir
        )
}
