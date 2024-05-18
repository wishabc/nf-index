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
    bedmap --sweep-all \
            --delim "\t" \
            --max ${params.index_file} tmp.bg \
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

    paste ${columns} > matrix.txt
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