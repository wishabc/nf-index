process signal_proportion_in_peaks {
    tag "${ag_id}"
	//label 'high_mem'
	publishDir "${params.outdir}/signal_prop/"

    input:
        tuple val(ag_id), path(peaks_file), path(cutcounts_file)
    
    output:
        tuple val(ag_id), path(name)

    script:
    name = "${ag_id}.prop_signal_stats.txt"
    """
    echo -e "length\tcutcounts\tdistance\tsample_id" > ${name}

    IFS=',' read -r -a distance_array <<< "${params.distances}"

    for DISTANCE in "\${distance_array[@]}"; do
        bash $moduleDir/bin/signal_prop_at_distance.sh \
            ${ag_id} \
            \$DISTANCE \
            ${peaks_file} \
            ${cutcounts_file} \
            ${params.chrom_sizes} >> ${name}
    done
    """
}


workflow {
    params.distances = "0,10,50,150,500,1000,5000,10000"
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.index_peaks_file),
            file(row.cutcounts_file),
        ))
        | signal_proportion_in_peaks
        | map(it -> it[1])
        | collectFile(
            storeDir: params.outdir,
            skip: 1,
            keepHeader: true,
            name: "signal_proportion_at_distance.txt"
        ) 
}