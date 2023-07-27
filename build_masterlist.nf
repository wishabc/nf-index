


process collate_and_chunk {
    conda params.conda
    
    input:
        path peak_files
    
    output:
        path "chunk*.bed"

    script:
    """
    
    for f in ${peak_files}; 
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
    module "R/3.3.3"

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
        ${chunk_file.parent} 
    """
}


process resolve_overlaps {

    module "R/3.3.3"

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
    module "kentutil:bedops"
    publishDir params.outdir

    input:
        path "DHSs_all/*"
        path "DHSs_nonovl_core/*"
        path "DHSs_nonovl_any/*"
    
    output:
        path "${prefix}*"
    
    script:
    prefix = "masterlist"
    """
    bash $moduleDir/bin/code_gen_masterlist.sh \
        ${prefix} \
        ${params.chrom_sizes_bed} \
        ./
    """
}

workflow {
    chunks = Channel.fromPath(params.peaks_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> it.peaks)
        | collate_and_chunk
        | process_chunk
        | resolve_overlaps
    
    merge_chunks(
        chunks[0].collect(sort: true), 
        chunks[1].collect(sort: true), 
        chunks[2].collect(sort: true)
    )
}