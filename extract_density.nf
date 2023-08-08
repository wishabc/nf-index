

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


workflow {	
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(it -> tuple(it.ag_id, file(it.normalized_density_file)))
        | extract_max_density
}