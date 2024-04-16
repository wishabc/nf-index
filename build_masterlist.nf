#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process collate_and_chunk {
    conda params.conda
    label "highmem"
    
    input:
        path "peaks/peaks*.starch"
    
    output:
        path "chunk*.bed"

    script:
    """
    
    for f in peaks/*; 
        do unstarch \$f >> tmp.bed;
    done

    sort-bed tmp.bed \
        | awk -v minChunkSep=10000 \
            -v minChunkSep2=4000 \
            -v minPerChunk=500 \
            -v maxPerChunk=100000 \
            -f $moduleDir/bin/index_scripts/chunk_bed.awk
    """
}


process process_chunk {
    conda params.conda
    tag "${prefix}"
    label "medmem"

    input:
        path chunk_file
    
    output:
        path "DHSs_all/${prefix}.bed"
        path "peaks_all/${prefix}.bed"
    
    script:
    prefix = "${chunk_file.baseName}"
    """
    Rscript $moduleDir/bin/index_scripts/code_build.R \
        ${prefix} \
        $moduleDir/bin/index_scripts \
        ./
    """
}


process resolve_overlaps {
    conda params.conda
    tag "${prefix}"
    label "medmem"

    input:
        path chunk_file, name: "DHSs_all/*"
        path peaks_file, name: "peaks_all/*"
    
    output:
        path "DHSs_all/${prefix}.bed"
        path "DHSs_nonovl_core/${prefix}.bed"
        path "DHSs_nonovl_any/${prefix}.bed"
    
    script:
    prefix = "${chunk_file.baseName}"
    """
    Rscript $moduleDir/bin/index_scripts/code_overlap.R \
        ${prefix} \
        $moduleDir/bin/index_scripts
    """
}

process merge_chunks {
    conda params.conda
    publishDir "${params.outdir}/unfiltered_masterlists"
    label "highmem"
    scratch true

    input:
        path filepaths_all 
        path filepaths_nonovl_core
        path filepaths_nonovl_any
    
    output:
        path "masterlist*${params.masterlist_id}*", emit: all
        path "masterlist_DHSs_${params.masterlist_id}_all_chunkIDs.bed", emit: non_merged
    
    script:
    """
    # workaround for too big .command.run
    mkdir "DHSs_all/"
    while read line; do
        ln -s \$line \$PWD/DHSs_all/
    done < ${filepaths_all}

    mkdir "DHSs_nonovl_core/"
    while read line; do
        ln -s \$line \$PWD/DHSs_nonovl_core/
    done < ${filepaths_nonovl_core}

    mkdir "DHSs_nonovl_any/"
    while read line; do
        ln -s \$line \$PWD/DHSs_nonovl_any/
    done < ${filepaths_nonovl_any}

    cat ${params.chrom_sizes} \
        | awk -v OFS='\t' '{ print \$1,0,\$2 }' \
        > chrom_sizes.bed
    
    bash $moduleDir/bin/index_scripts/code_gen_masterlist.sh \
        ${params.masterlist_id} \
        chrom_sizes.bed \
        ./
    """
}

workflow buildIndex {
    take:
        peaks
    main:
        chunks = peaks
            | collect(sort: true)
            | collate_and_chunk
            | flatten()
            | process_chunk
            | resolve_overlaps
 
        index = merge_chunks(
            // workaround
            chunks[0].map(it -> it.toString()).collectFile(name: 'all.paths.txt', newLine: true), 
            chunks[1].map(it -> it.toString()).collectFile(name: 'no_core.paths.txt', newLine: true), 
            chunks[2].map(it -> it.toString()).collectFile(name: 'no_any.paths.txt', newLine: true)
        ).non_merged
    emit:
        index
        process_chunk.out[1]
}

workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> file(row.peaks_file))
        | buildIndex
}
