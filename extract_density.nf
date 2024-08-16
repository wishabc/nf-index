#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process extract_max_density {
    conda params.conda
    publishDir "${params.outdir}/density"
    tag "${ag_id}"
    scratch true

    input:
        tuple val(ag_id), path(peaks_file)
    
    output:
        tuple val(ag_id), path(density)
    
    script:
    density = "${ag_id}.mean_max.tsv"
    """
    bigWigToBedGraph ${peaks_file} tmp.bg 
    cat tmp.bg \
        | awk -v OFS='\t' '{print \$1,\$2,\$3,"${ag_id}",\$4}' \
        | bedmap --sweep-all \
            --delim "\t" \
            --max ${params.index_file} - \
        > ${density}
    """
}

process collect_matrix {
    conda params.conda
    publishDir params.outdir
    label "bigmem"
    scratch true

    input:
        path columns
    
    output:
        path matrix
        path "samples_order.txt"
    
    script:
    matrix = "matrix.density.npy"
    """
    echo "${columns}" \
        | tr " " "\n"  \
		| xargs -I file basename file \
		| cut -d. -f1 \
		| tr "\n" "\t" > samples_order.txt

    paste ${columns} | sed 's/\\<NAN\\>/0/g' > matrix.txt
    python3 $moduleDir/bin/convert_to_numpy.py matrix.txt ${matrix}
    """

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

workflow {	
    matrix = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ag_id, file(row.normalized_density_bw)))
        | extract_max_density
        | map(it -> it[1])
        | collect(sort: true)
        | collect_matrix
    
}


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


process create_genome_chunks {
	memory 500.MB
	conda "${params.conda}"
	scratch true

	output:
		path "genome_chunks.bed"

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
            }' > genome_chunks.bed
	"""
}

process apply_wiggletools {
    conda params.conda
    publishDir "${params.outdir}"
    scratch true

    input:
        tuple val(function), val(chunk)
        path bigwigs
    
    output:
        tuple val(function), path(name)
    
    script:
    name = "normalized.${function}.${chunk}.tsv"
    """
    wiggletools write normalized.${function}.wig \
        ${function} \
        seek ${chunk} \
        `cat ${bigwigs}` > ${name}
    """
}


workflow averageTracks {
    functions = Channel.of('median', 'mean', 'max')
    functions_and_chunks = create_genome_chunks()
        | flatMap(n -> n.split())
        | map(it -> it.replaceAll(':', ' ').replaceAll('-', ' '))
        | combine(functions) // chunk, function

    bigwigs = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> file(row.normalized_density_bw))
        | collectFile(sort: true, name: 'bigwigs.txt', newLine: true)
     
    apply_wiggletools(functions_and_chunks, bigwigs)
        | collectFile(
                storeDir: params.outdir,
                sort: true,
            ) { it -> [ "normalized.${it[0]}.tsv", it[1].text ] }
}