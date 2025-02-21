#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process variance_partition {

    tag "${start_index}"
    scratch true
    label "medmem"
    conda params.conda

    input:
        tuple val(start_index), path(norm_matrix), path(masterlist), path(samples_meta), val(formula)
    
    output:
        path name
    
    script:
    name = "${start_index}.variance_partition.tsv"
    """
    Rscript $moduleDir/bin/post_processing/variance_partition.R \
        ${samples_meta} \
        ${start_index} \
        ${params.chunk_size} \
        ${norm_matrix} \
        ${masterlist} \
        '${formula}' \
        ${name} 
    """
}

process sort_bed {

    conda params.conda
    publishDir params.outdir
    label "medmem"

    input:
        path unsorted_bed
    
    output:
        path name
    
    script:
    name = "masterlist.vp_annotated.sorted.bed"
    """
    head -1 ${unsorted_bed} > ${name}
    tail -n+2 ${unsorted_bed} | sort-bed - >> ${name}
    """

}

workflow variancePartition {
    take:
        data // normalized_matrix, masterlist, samples_file, formula
    main:
        out = data
            | flatMap(it -> (1..it[1].countLines()))
            | collate(params.chunk_size, remainder=true)
            | map(it -> it[0])
            | combine(data) // chunk_start, normalized_matrix, masterlist, samples_order, samples_file, formula
            | variance_partition
            | collectFile(
                name: "masterlist.vp_annotated.bed",
                keepHeader: true,
                sort: true,
                skip: 1
            )
            | sort_bed
            | combine(data.map(it -> it[3]))
    emit:
        out  
}


// workflow {
//     params.h5file = "${params.outdir}/matrices.h5"
//     variancePartition(
//         Channel.fromPath(params.masterlist),
//         Channel.fromPath(params.h5file)
//     )
// }
