#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process collate_and_chunk {
    conda params.conda
    
    input:
        path "peaks/peaks*.ext"
    
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
            -f $moduleDir/bin/chunk_bed.awk
    """
}


process process_chunk {
    conda params.conda
    tag "${prefix}"

    input:
        path chunk_file
    
    output:
        path "DHSs_all/${prefix}.bed"
        path "peaks_all/${prefix}.bed"
    
    script:
    prefix = "${chunk_file.baseName}"
    """
    Rscript $moduleDir/bin/code_build.R \
        ${prefix} \
        $moduleDir/bin \
        ./
    """
}


process resolve_overlaps {
    conda params.conda
    tag "${prefix}"

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
    Rscript $moduleDir/bin/code_overlap.R \
        ${prefix} \
        $moduleDir/bin
    """
}

process merge_chunks {
    conda params.conda
    publishDir "${params.outdir}/unfiltered_masterlists"
    scratch true

    input:
        path "DHSs_all/*"
        path "DHSs_nonovl_core/*"
        path "DHSs_nonovl_any/*"
    
    output:
        path "masterlist*${params.masterlist_id}*", emit: all
        path "masterlist_DHSs_${params.masterlist_id}_all_chunkIDs.bed", emit: non_merged
    
    script:
    """
    cat ${params.chrom_sizes} \
        | awk -v OFS='\t' '{ print \$1,0,\$2 }' \
        > chrom_sizes.bed
    
    bash $moduleDir/bin/code_gen_masterlist.sh \
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
 
        masterlist = merge_chunks(
            chunks[0].collect(sort: true), 
            chunks[1].collect(sort: true), 
            chunks[2].collect(sort: true)
        ).non_merged
            | filter_masterlist // returns tuple(masterlist, mask)
            | map(it -> it[0]) // masterlist
    emit:
        masterlist
	
}

workflow {
    Channel.fromPath(params.peaks_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> it.peaks)
        | buildIndex
}
