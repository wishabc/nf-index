
process convert_to_numpy {
    conda params.conda
    publishDir "${params.outdir}/raw_matrices"
    label "highmem"

    input:
        tuple val(prefix), path(matrix)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.raw.matrix.npy"
    dtype = prefix.contains('binary') ? 'bool' : (prefix.contains('counts') ? 'int' : 'float')
    """
    python3 $moduleDir/bin/convert_to_numpy.py \
        ${matrix} \
        ${name} \
        --dtype ${dtype}
    """
}

process apply_filter {
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
	"""
    python3 $moduleDir/bin/apply_mask.py \
        ${matrix} \
        ${name} \
        ${mask}
	"""
}


process filter_masterlist {
    conda params.conda

    // publishDir "${params.outdir}", pattern: "${only_autosomes_masterlist}"
    // publishDir "${params.outdir}/masks/masterlist_no_header", pattern: "${filtered_masterlist}"
    publishDir "${params.outdir}/masks", pattern: "*.mask.txt"

    input:
        tuple val(prefix), path(binary_matrix), path(masterlist, name: 'masterlist.bed')
    
    output:
	    tuple path(non_zero_rows), path(filtered_mask), path(only_autosomes_mask)

    script:
    non_zero_rows = "${prefix}.non_zero_rows.mask.txt"
    filtered_masterlist = "${prefix}_DHSs.blacklistfiltered.bed"
    filtered_mask = "${prefix}.filtered_DHS.mask.txt"
    only_autosomes_mask = "${prefix}.filtered.autosomes.mask.txt"
    only_autosomes_masterlist = "${prefix}.only_autosomes.filtered.bed"
    """
    zcat ${binary_matrix} \
        | awk '{ if (/1/) print 1; else print 0; }' > ${non_zero_rows}
    
    bedmap --bases ${masterlist} \
        ${params.encode_blacklist_regions} \
        |  awk -F'\t' '{ if(\$1 > 0) print 1; else print 0}' \
        > blacklist_rows.txt

    python3 $moduleDir/bin/filter_dhs.py \
        ${masterlist} \
        ${non_zero_rows} \
        blacklist_rows.txt \
        ${filtered_masterlist} \
        ${filtered_mask} \
        --singletons_strategy ${params.singletons_strategy}

    cat ${masterlist} \
		| awk '{print (\$1 ~ /^chr[0-9]+/) ? 1 : 0}' \
		> autosomes.mask.txt
    
    awk 'NR==FNR {mask1[NR]=\$0; next} \
        {print mask1[FNR] * \$0}' \
        ${filtered_mask} \
        autosomes.mask.txt > ${only_autosomes_mask}

	awk 'NR==FNR {mask[NR]=\$0; next} mask[FNR] == 1' \
		${only_autosomes_mask} ${masterlist} > ${only_autosomes_masterlist}
    """
}

workflow filterAndConvertToNumpy {
    take:
        masterlist
        raw_matrices

    main:
        masterlist_and_mask = raw_matrices 
            | filter { it[0].startsWith('binary') }
            | combine(masterlist)
            | filter_masterlist

        raw_np = raw_matrices
            | convert_to_numpy

        // w_autosomes_matrices = raw_np 
        //     | combine(
        //         masterlist_and_mask.map(it -> it[2])
        //     )

        // npy_matrices = raw_np
        //     | map(it -> tuple("${it[0]}.only_autosomes", it[1]))
        //     | combine(
        //         masterlist_and_mask.map(it -> it[4])
        //     )
        //     | mix(w_autosomes_matrices)
        //     | apply_filter
    emit:
        masterlist_and_mask
        raw_np
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

workflow getMasks {
    params.base_dir = params.outdir

    matrices = Channel.fromPath("${params.base_dir}/raw_matrices/matrix.${it}.mtx.gz")
        | combine(params.index_file)
        | filter_masterlist
}
