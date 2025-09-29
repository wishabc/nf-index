#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process variance_partition {

    tag "${start_index}"
    scratch true
    label "medmem"
    conda params.conda

    input:
        tuple val(start_index), val(prefix), path(norm_matrix), path(samples_order), path(masterlist), path(samples_meta), val(variance_partition_formula)
    
    output:
        tuple val(prefix), path(name)
    
    script:
    name = "${prefix}.${start_index}.variance_partition.tsv"
    """
    Rscript $moduleDir/bin/r_scripts/variance_partition.R \
        ${samples_order} \
        ${samples_meta} \
        ${start_index} \
        ${params.chunk_size} \
        ${norm_matrix} \
        ${masterlist} \
        '${variance_partition_formula}' \
        ${name} 
    """
}

process sort_bed {

    conda params.conda
    publishDir params.outdir
    label "medmem"

    input:
        tuple val(prefix), path(unsorted_bed)
    
    output:
         tuple val(prefix), path(name)
    
    script:
    name = "${prefix}.vp_annotated.sorted.bed"
    """
    head -1 ${unsorted_bed} > ${name}
    tail -n+2 ${unsorted_bed} | sort-bed - >> ${name}
    """

}

workflow variancePartition {
    take:
        data // prefix, vst_matrix, samples_order, masterlist, samples_file
    main:
        out = data
            | flatMap(it -> (1..it[3].countLines()))
            | collate(params.chunk_size, remainder=true)
            | map(it -> it[1]) // chunk_start
            | combine(data) // chunk_start, prefix, vst_matrix, samples_order, masterlist, samples_file, variance_partition_formula
            | variance_partition // prefix, vp_annotated_chunk
            | collectFile(
                keepHeader: true,
                sort: true,
                skip: 1
            ) { [
                "${it[0]}.masterlist.vp_annotated.bed", // name
                it[1] // content
            ] }
            | map(it -> tuple(it[0].name.replaceAll(".masterlist.vp_annotated.bed", ""), it[1])) // prefix, vp_annotated_masterlist
            | sort_bed
    emit:
        out  
}


workflow {
    input_data = Channel.of(
        tuple(
            file(params.norm_matrix),
            file(params.masterlist),
            file(params.samples_meta),
            params.variance_partition_formula
        )
    ) | variancePartition
}
