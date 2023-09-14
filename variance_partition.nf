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
    Rscript $moduleDir/bin/variance_partition.R \
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
    publishDir params.oudir
    label "highmem"
    
    input:
        path binary_matrix
        path vst_matrix
        path samples_file

    output:
        path name
    script:
    name = "matrices.h5"
    """
    python3 $moduleDir/bin/convert_to_h5.py \
        ${vst_matrix} \
        ${samples_names} \
        ${name} \
        --binary ${binary_matrix}
    """
}

workflow variancePartition {
    take:
        masterlist
        h5file
    main:
        out = masterlist
            | flatMap(it -> (1..it.countLines()))
            | collate(params.chunk_size)
            | map(it -> it[0])
            | combine(masterlist)
            | combine(h5file)
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
    params.chunk_size = 5000
    params.h5file = "$launchDir/${params.outdir}/matrices.h5"
    params.masterlist = "$launchDir/${params.outdir}/masterlist.filtered.bed"
    params.formula = "~ (1 | extended_annotation) + (1 | ln_finished_date) + (1 | frac_method) + (1 | is_primary_tissues_all)"
    variancePartition(
        Channel.fromPath(params.masterlist),
        Channel.fromPath(params.h5file)
    )
}
