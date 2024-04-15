
process apply_filter_and_convert_to_np {
	publishDir "${publishDirectory}", pattern: "${name}"
	label "highmem"
    tag "${prefix}"
	conda params.conda

	input:
		tuple val(prefix), path(matrix), path(mask)
	
	output:
		tuple val(prefix), path(name)
	
	script:
	name = "${prefix}.filtered.matrix.npy"
    publishDirectory = prefix.contains("only_autosomes") ?  "${params.outdir}" : "${params.outdir}/annotations"
    dtype = prefix.contains('binary') ? 'bool' : (prefix.contains('signal') ? 'int' : 'float')
	"""
    python3 $moduleDir/bin/convert_to_numpy.py \
        ${matrix} \
        ${name} \
        --dtype ${dtype} \
        --mask ${mask}
	"""
}


process filter_masterlist {
    conda params.conda

    publishDir "${params.outdir}", pattern: "${filtered_masterlist}"
    publishDir "${params.outdir}/masks", pattern: "*.mask.txt"

    input:
        tuple path(binary_matrix), path(masterlist)
    
    output:
	    tuple path(name), path(mask), path(filtered_masterlist), path(autosomes_mask), path(non_zero_rows)

    script:
    prefix = "masterlist"
    non_zero_rows = "${prefix}.non_zero_rows.mask.txt"
    name = "${prefix}_DHSs.blacklistfiltered.bed"
    mask = "${prefix}.bad_dhs.mask.txt"
    autosomes_mask = "${prefix}.filtered.autosomes.mask.txt"
    filtered_masterlist = "${prefix}.only_autosomes.filtered.bed"
    """
    zcat ${binary_matrix} | awk '{ if (/1/) print 1; else print 0; }' > ${non_zero_rows}
    bedmap --bases ${masterlist} \
        ${params.encode_blacklist_regions} \
        |  awk -F'\t' '{ if(\$1 > 0) print 1; else print 0}' \
        > blacklist_rows.txt

    python3 $moduleDir/bin/filter_dhs.py \
        ${prefix} \
        .5 \
        blacklist_rows.txt \
        ${non_zero_rows} \
        ${masterlist} \
        ${name} \
        ${mask}
    
    cat ${masterlist} \
		| awk '{print (\$1 ~ /^chr[0-9]+/) ? 1 : 0}' \
		> autosomes.mask.txt
    
    awk 'NR==FNR {mask1[NR]=\$0; next} \
        {print mask1[FNR] * \$0}' ${mask} autosomes.mask.txt > ${autosomes_mask}

	awk 'NR==FNR {mask[NR]=\$0; next} mask[FNR] == 1' \
		${autosomes_mask} ${masterlist} > ${filtered_masterlist}
    """
}

workflow filterAndConvertToNumpy {
    take:
        masterlist
        raw_matrices
    main:
        masterlist_and_mask = raw_matrices 
            | filter { it[0] == 'binary'}
            | map(it -> it[1])
            | combine(masterlist)
            | filter_masterlist

        w_autosomes_matrices = raw_matrices 
            | combine(masterlist_and_mask.map(it -> it[1]))

        npy_matrices = raw_matrices
            | map(it -> tuple("${it[0]}.only_autosomes", it[1]))
            | combine(
                masterlist_and_mask.map(it -> it[3])
            )
            | mix(w_autosomes_matrices)
            | apply_filter_and_convert_to_np
    emit:
        masterlist_and_mask
        npy_matrices
}

workflow {
    params.base_dir = params.outdir

    matrices = Channel.of('binary', 'counts')
        | map(it -> tuple(it, file("${params.base_dir}/raw_matrices/matrix.${it}.mtx.gz")))

    filterAndConvertToNumpy(
        Channel.fromPath(params.index_file),
        matrices
    )
}
