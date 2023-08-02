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
    // FIXME remove module
    module "kentutil:bedops"
    publishDir "${params.outdir}/unfiltered_masterlists"

    input:
        path "DHSs_all/*"
        path "DHSs_nonovl_core/*"
        path "DHSs_nonovl_any/*"
    
    output:
        path "masterlist*${params.masterlist_id}*", emit: all
        path "masterlist_DHSs_${params.masterlist_id}_all_chunkIDs.bed", emit: non_merged
    
    script:
    """
    bash $moduleDir/bin/code_gen_masterlist.sh \
        ${params.masterlist_id} \
        ${params.chrom_sizes_bed} \
        ./
    """
}

process filter_masterlist {
    publishDir params.outdir
    conda params.conda

    input:
        path masterlist
    
    output:
	    path name

    script:
    prefix = "masterlist"
    name = "masterlist_DHSs.blacklistfiltered.bed"
    """
    bedmap --bases ${masterlist} ${params.encode_blacklist_regions} \
        |  awk -F'\t' '{ if(\$1 > 0) print (NR-1)}' \
        > blacklist_rows.txt

    python3 $moduleDir/bin/DHS_filter.py \
        ${prefix} \
        .5 \
        blacklist_rows.txt \
        ${masterlist} \
        ${name}
    """
}

process annotate_masterlist {
    publishDir params.outdir

    input: 
        path filtered_masterlist

    output:
        path name

    script:
    name = "masterlist_DHSs_${params.masterlist_id}.filtered.annotated.bed"
    """
    bash $moduleDir/bin/simpleAnnotations.sh \
        ${filtered_masterlist} \
        ${params.encode3} \
        ${params.gencode} \
        ${params.gwas_catalog}
    
    echo -e "#chr\tstart\tend\tdhs_id\ttotal_signal\tnum_samples\tnum_peaks\tdhs_width\tdhs_summit\tcore_start\tcore_end\tmean_signal" > masterlist_header.txt
    echo -e "is_encode3\tencode3_ovr-fraction\tdist_tss\tgene_name\tnum_gwasCatalog_variants" > simpleAnnotations_header.txt

    paste masterlist_header.txt simpleAnnotations_header.txt > header.txt
    paste ${filtered_masterlist} \
        is_encode3.txt \
        dist_gene.txt \
        gwas_catalog_count.txt \
        | cat header.txt - > ${name}
 
    """
}


workflow {
    chunks = Channel.fromPath(params.peaks_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> it.peaks)
        | collect(sort: true)
        | collate_and_chunk
        | flatten()
        | process_chunk
        | resolve_overlaps
    
    masterlists = merge_chunks(
        chunks[0].collect(sort: true), 
        chunks[1].collect(sort: true), 
        chunks[2].collect(sort: true)
    ).non_merged
        | filter_masterlist
        | annotate_masterlist
	
    	
}
