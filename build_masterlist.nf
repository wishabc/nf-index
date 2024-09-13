#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { filterAndConvertToNumpy } from "./filter_peaks"
include { convert_index_to_anndata } from "./convert_to_anndata"

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

process get_samples_order {

    publishDir params.outdir
    
    output:
        path name

    script:
    name = "samples_order.txt"
    """
    awk -F"\t" -v col="ag_id" \
        'NR==1{for(i=1;i<=NF;i++)if(\$i==col)c=i}NR>1{if(c)print \$c}' \
            ${params.samples_file} > ${name}
    """
}


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
    scratch true

    input:
        tuple path(chunk_file), path(chunks_order), path(samples_order)

    output:
        tuple val(prefix), path(binary)

    script:
    prefix = "binary.index"
    binary = "${chunk_file.baseName}.${prefix}.txt"
    """
    $moduleDir/bin/index_scripts/writeRowsPerChunkForMatrices \
        ${samples_order} \
        ${chunks_order} \
        ${chunk_file} \
        tmp.txt \
        ${binary}
    """
}


process collect_chunks {
    publishDir "${params.outdir}/raw_matrices"
    tag "${prefix}"

    input:
        tuple val(prefix), path("chunks/*")
    
    output:
        tuple val(prefix), path(matrix)

    script:
    matrix = "${prefix}.raw.matrix.mtx.gz"
    """
    ls chunks/ \
        | wc -l \
        | awk -v dir="chunks/" \
            '{ \
                for(i=1;i<=\$1;i++){ \
                    printf("%s/chunk%04d.${prefix}.txt ",dir,i) \
                } printf("\\n"); \
            }' \
        | xargs cat \
        | gzip > ${matrix}
    """
}

process annotate_masterlist {
    conda params.conda
    publishDir "${params.outdir}", pattern: "${name}"
    scratch true
    label "highmem"
    errorStrategy 'ignore'

    input: 
        tuple path(binary_matrix), path(samples_order), path(masterlist)

    output:
        tuple path(binary_matrix), path(samples_order), path(name)

    script:
    name = "masterlist_DHSs_all_chunks.${params.masterlist_id}.annotated.bed"
    """
     python $moduleDir/bin/annotations/spot1Annotations.py \
        ${binary_matrix} \
        ${samples_order} \
        ${params.samples_file} \
        spot1_metrics.txt
 
    bash $moduleDir/bin/annotations/simpleAnnotations.sh \
        ${masterlist} \
        ${params.encode3} \
        ${params.gwas_catalog} \
	    ${params.repeats} \
        simple_annotations.txt
    
    bash $moduleDir/bin/annotations/gencodeAnnotations.sh \
        ${masterlist} \
        ${params.gencode} \
        ${params.chrom_sizes} \
        gencode_annotations.txt

    bash $moduleDir/bin/annotations/gcContentAnnotations.sh \
        ${masterlist} \
        ${params.genome_fasta} \
        ${params.mappable_file} \
        gc_content.txt

    echo -e "#chr\tstart\tend\tdhs_id\ttotal_signal\tnum_samples\tnum_peaks\tdhs_width\tdhs_summit\tcore_start\tcore_end" \
        | cat - ${masterlist} \
        | paste - \
            simple_annotations.txt \
            gencode_annotations.txt \
            spot1_metrics.txt \
            gc_content.txt > ${name}
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
            chunks[0].map(it -> it.toString()).collectFile(name: 'all.paths.txt', newLine: true, sort: true), 
            chunks[1].map(it -> it.toString()).collectFile(name: 'no_core.paths.txt', newLine: true, sort: true), 
            chunks[2].map(it -> it.toString()).collectFile(name: 'no_any.paths.txt', newLine: true, sort: true)
        ).non_merged

        chunks_order = index
            | get_chunks_order

        samples_order = get_samples_order()

        binary_matrix = process_chunk.out[1]
            | combine(chunks_order)
            | combine(samples_order)
            | write_rows
            | groupTuple()
            | collect_chunks // prefix, matrix
            | combine(index)
            | filterAndConvertToNumpy
        
        out = binary_matrix[1] // prefix, np_matrix
            | map(it -> it[1])
            | combine(samples_order)
            | combine(index)
            

    emit:
        out
        binary_matrix[0]
}

workflow annotateMasterlist {
    Channel.of(
        tuple(
            file("${params.index_dir}/raw/binary.index.raw.matrix.npy"),
            file("${params.index_dir}/samples_order.txt"),
            file("${params.index_dir}/unfiltered_masterlist/masterlist_DHSs_${params.masterlist_id}_all_chunkIDs.bed")
        )
    ) | annotate_masterlist
}

workflow {
    index_data = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> file(row.peaks_file))
        | buildIndex
    
    masks = index_data[1] | collect(sort: true, flat: true)

    annotated_masterlist = index_data[0]
        | annotate_masterlist // matrix, samples_order, annotated_index
        | combine(masks)
        | convert_index_to_anndata
}
