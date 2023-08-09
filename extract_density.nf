

process extract_max_density {
    conda params.conda
    publishDir "${params.outdir}/density"
    tag "${ag_id}"

    input:
        tuple val(ag_id), path(peaks_file)
    
    output:
        tuple val(ag_id), path(density)
    
    script:
    density = "${ag_id}.mean_max.tsv"
    """
    bedmap --sweep-all \
        --delim "\t" \
        --max ${params.index_file} ${peaks_file} \
    > ${density}
    """
}

process collect_matrix {
    conda params.conda
    publishDir params.outdir
    label "bigmem"

    input:
        path columns
    
    output:
        path matrix
        path "samples_order.txt"
    
    script:
    matrix = "matrix.density.npy"
    """
    echo "${columns}" | tr " " "\n"  \
		| xargs -I file basename file \
		| cut -d. -f1 \
		| tr "\n" "\t" > samples_order.txt

    paste ${columns} > matrix.txt
    python3 $moduleDir/bin/convert_to_numpy.py matrix.txt ${matrix}
    """

}


workflow {	
    matrix = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> tuple(it.ag_id, file(it.normalized_density_file)))
        | extract_max_density
        | map(it -> it[1])
        | collect(sort: true)
        | collect_matrix
    
}