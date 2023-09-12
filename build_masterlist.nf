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
    
    bash $moduleDir/bin/code_gen_masterlist.sh \
        ${params.masterlist_id} \
        chrom_sizes.bed \
        ./
    """
}


process filter_masterlist {
    conda params.conda

    input:
        path masterlist
    
    output:
	    tuple path(name), path(mask)

    script:
    prefix = "masterlist"
    name = "${prefix}_DHSs.blacklistfiltered.bed"
    mask = "${prefix}.mask.txt"
    """
    bedmap --bases ${masterlist} ${params.encode_blacklist_regions} \
        |  awk -F'\t' '{ if(\$1 > 0) print (NR-1)}' \
        > blacklist_rows.txt

    python3 $moduleDir/bin/DHS_filter.py \
        ${prefix} \
        .5 \
        blacklist_rows.txt \
        ${masterlist} \
        ${name} \
        ${mask}
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
            // workaround
            chunks[0].map(it -> it.toString()).collectFile(name: 'all.paths.txt', newLine: true), 
            chunks[1].map(it -> it.toString()).collectFile(name: 'no_core.paths.txt', newLine: true), 
            chunks[2].map(it -> it.toString()).collectFile(name: 'no_any.paths.txt', newLine: true)
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


// DEFUNC
process get_chunks_order {
    input:
        path masterlist
    
    output:
        path name

    script:
    name = "chunks_order.txt"
    """
    cut -f4 ${masterlist} > ${name}
    """
}
process write_rows {

    tag "${chunk_file.simpleName}"
    input:
        tuple path(chunk_file), path(chunks_order), path(samples_order)

    output:
        path signal, emit: signal
        path binary, emit: binary

    script:
    signal = "${chunk_file.baseName}.signal.txt"
    binary = "${chunk_file.baseName}.binary.txt"
    """
    /net/seq/data/projects/SuperIndex/erynes/masterLists/writeRowsPerChunkForMatrices \
        ${samples_order} \
        ${chunks_order} \
        ${chunk_file} \
        ${signal} \
        ${binary}
    """
}

process collect_chunks {

    publishDir params.outdir
    tag "${prefix}"

    input:
        tuple val(prefix), path("chunks/*")
    
    output:
        path matrix

    script:
    matrix = "matrix.${prefix}.mtx.gz"
    """
    ls chunks/ | wc -l \
        | awk -v dir="chunks/" \
         '{for(i=1;i<=\$1;i++){printf("%s/chunk%04d.${prefix}.txt ",dir,i)}printf("\\n");}' \
        | xargs cat | gzip > ${matrix}
    """
}

workflow createMatrices {
    chunks_order = Channel.fromPath("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/output/masterlist_DHSs_0802.filtered.annotated.bed")
        | get_chunks_order

    samples_order = Channel.of(file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/${params.outdir}/samples_order.txt"))
    rows = Channel.fromPath("/net/seq/data2/projects/sabramov/SuperIndex/dnase-wouter-style-matrices/peaks_list.txt")
        | splitCsv(header: false)
        | map(it -> it[0])
        | combine(chunks_order)
        | combine(samples_order)
        | write_rows
    
    binary_data = rows.binary
        | collect(sort: true)
        | map(it -> tuple("binary", it))
    
    rows.signal
        | collect(sort: true)
        | map(it -> tuple("signal", it))
        | mix(binary_data)
        | collect_chunks
    // file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/work/tmp/2a/10c57903166faeb052d56a4ace1a68/no_core.paths.txt")
    // file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/work/tmp/93/23361e11b94264bc72ccbac2af1af4/no_any.paths.txt")
}