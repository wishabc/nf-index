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
        ${params.formula} \
        ${name} 
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
    params.filtered_masterlist = "$launchDir/${params.outdir}/masterlist.filtered.bed"
    params.formula = "~ dedupped_subsampled_spot1 + log(read_depth) + dupRate_5M + (1 | donor_sex) + (1 | library_kit) + (1 | short_ontology)"
    variancePartition(
        Channel.fromPath(params.filtered_masterlist),
        Channel.fromPath(params.h5file)
    )
}
