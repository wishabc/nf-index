#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process variance_partition {
    conda "/home/sabramov/miniconda3/envs/condR-clone"
    tag "${start_index}"
    scratch true

    input:
        tuple val(start_index), path(masterlist), path(h5file)
    
    output:
        path name
    
    script:
    name = "${start_index}.variance_partition.tsv"
    end_index = start_index + params.chunk_size - 1
    """
    Rscript $moduleDir/bin/variance_partition/variance_partition.R \
        ${params.samples_file} \
        ${start_index} \
        ${params.chunk_size} \
        ${h5file} \
        ${masterlist} \
        '${params.formula}' \
        ${name} 
    """
}

process convert_to_h5 {
    conda params.conda
    publishDir params.outdir
    label "highmem"
    
    input:
        path binary_matrix
        path vst_matrix
        path samples_names

    output:
        path name
    script:
    name = "matrices.h5"
    """
    python3 $moduleDir/bin/variance_partition/convert_to_h5.py \
        ${vst_matrix} \
        ${samples_names} \
        ${name} \
        --binary ${binary_matrix}
    """
}
process sort_bed {

    conda params.conda
    publishDir params.outdir

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
        masterlist
        h5file
    main:
        params.chunk_size = 5000
        out = masterlist
            | flatMap(it -> (1..it.countLines()))
            | collate(params.chunk_size)
            | map(it -> it[0])
            | combine(masterlist)
            | combine(h5file)
            | variance_partition
            | collectFile(
                name: "masterlist.vp_annotated.bed",
                keepHeader: true,
                sort: true,
                skip: 1
            )
            | sort_bed
    emit:
        out  
}

params.masterlist = "${params.outdir}/masterlist.only_autosomes.filtered.bed"

workflow convertToH5 {
    h5 = convert_to_h5(
        Channel.fromPath("${params.outdir}/binary.only_autosomes.filtered.matrix.npy"),
        Channel.fromPath("${params.outdir}/deseq_normalized.only_autosomes.filtered.sf.vst.npy"),
        Channel.fromPath("${params.outdir}/samples_order.txt")
    )
    variancePartition(Channel.fromPath(params.masterlist), h5)
}

workflow {
    params.h5file = "${params.outdir}/matrices.h5"
    params.formula = "~ (1 | core_annotation) + dup_rate + SPOT1_score + I(SPOT1_score^2)"
    variancePartition(
        Channel.fromPath(params.masterlist),
        Channel.fromPath(params.h5file)
    )
}
