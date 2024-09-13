#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process create_genome_chunks {
	memory 500.MB
	conda "${params.conda}"
	scratch true

	output:
		stdout

	script:
	"""
	cat ${params.chrom_sizes} \
	    | grep -v chrX \
        | grep -v chrY \
        | grep -v chrM \
        | grep -v _random \
        | grep -v _alt \
        | grep -v chrUn \
        | awk -v step=${params.chunksize} -v OFS="\t" \
            '{ \
                for(i=step; i<=\$2; i+=step) { \
                    print \$1"\t"i-step+1"\t"i; \
                } \
                print \$1"\t"i-step+1"\t"\$2; \
            }' \
        | sed 's/\t/_/g'
	"""
}

process collect_stats_for_chunk {
    conda params.conda
    publishDir "${params.outdir}/chunks"
    tag "${chunk_id}"
    errorStrategy 'retry'
    memory { 80.GB * task.attempt }


    input:
        val(chunk_id)
    
    output:
        path "*.${chunk_id}.bed"
    
    script:
    """
    python3 $moduleDir/bin/post_processing/stats_on_bw.py ${chunk_id} ${params.samples_file}
    """
}


workflow averageTracks {
    chunks = create_genome_chunks()
        | flatMap(n -> n.split())
        | collect_stats_for_chunk
        | flatten()
        | map(it -> tuple(it.simpleName, it))
        | collectFile(
                storeDir: params.outdir,
                sort: true,
            ) { it -> [ "normalized.${it[0]}.tsv", it[1].text ] }
}

// DEFUNC

process count_peaks {
    conda params.conda
    tag "${ag_id}"
    scratch true

    input:
        tuple val(ag_id), path(peaks_file)
    
    output:
        tuple val(ag_id), val(peaks_file.name), stdout
    
    script:
    """
    unstarch ${peaks_file} | wc -l
    """
}

workflow countPeaks {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.peaks_file)))
        | count_peaks
        | map(it -> "${it[0]}\t${it[1]}\t${it[2]}")
        | collectFile(
            name: "peaks_count.tsv",
            storeDir: params.outdir,
            newLine: true
        )
}


workflow tmp {	
    matrix = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.normalized_density_bw.replaceAll(".bw", ".starch"))))
        | extract_max_density
        | map(it -> it[1])
        | collect(sort: true)
        | collect_matrix
    
}