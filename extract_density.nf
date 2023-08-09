

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
    tail -n +2 ${params.index_file} \
        | bedmap --sweep-all \
            --delim "\t" \
             --max - ${peaks_file} \
        > ${density}
    """
}

process collect_matrix {
    conda params.conda
    publishDir params.outdir

    input:
        tuple val(ag_ids), path(columns)
    
    output:
        path matrix
    
    script:
    matrix = "matrix.density.tsv"
    """
    echo "${count_files}" | tr " " "\n"  \
		| xargs -I file basename file \
		| cut -d. -f1 \
		| tr "\n" "\t" > order.txt
	
	truncate -s -1 order.txt > samples_order.txt

    paste ${count_files} \
        | gzip -c > matrix.all.signal.txt.gz

    """

}


workflow {	
    matrix = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> tuple(it.ag_id, file(it.normalized_density_file)))
        | extract_max_density
        | collect(sort: true, flat: false)
        | collect_matrix
}