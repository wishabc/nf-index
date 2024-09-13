
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

    publishDir "${params.outdir}/masks"

    input:
        tuple val(prefix), path(binary_matrix), path(masterlist, name: 'masterlist.bed'), val(singletons_strategy)
    
    output:
	    tuple path(non_zero_rows), path(blacklist_rows), path(filtered_mask), path(only_autosomes_mask)

    script:
    non_zero_rows = "${prefix}.non_zero_rows.mask.txt"
    filtered_mask = "${prefix}.filtered_DHS.mask.txt"
    only_autosomes_mask = "${prefix}.filtered.autosomes.mask.txt"
    blacklist_rows = "${prefix}.blacklist_rows.mask.txt"
    """
    zcat ${binary_matrix} \
        | awk '{ if (/1/) print 1; else print 0; }' > ${non_zero_rows}
    
    grep -v "^#" ${masterlist} > masterlist.no_header.bed

    # FIXME to work with list of autosomes
    cat masterlist.no_header.bed \
		| awk '{print (\$1 ~ /^chr[0-9]+/) ? 1 : 0}' \
		> ${only_autosomes_mask}

    bedmap --bases masterlist.no_header.bed \
        ${params.encode_blacklist_regions} \
        |  awk -F'\t' \
            '{ if(\$1 > 0) print 1; else print 0}' \
        > ${blacklist_rows}

    python3 $moduleDir/bin/filter_dhs.py \
        masterlist.no_header.bed \
        ${non_zero_rows} \
        ${blacklist_rows} \
        ${filtered_mask} \
        --singletons_strategy ${singletons_strategy}
    """
}

workflow filterAndConvertToNumpy {
    take:
        raw_matrices

    main:
        masks = raw_matrices 
            | filter { it[0].startsWith('binary') }
            | filter_masterlist

        raw_np = raw_matrices
            | map(it -> tuple(it[0], it[1]))
            | convert_to_numpy

    emit:
        masks
        raw_np
}

workflow {
    params.base_dir = params.outdir
    singletons_strategy = Channel.of('keep_all')

    matrices = Channel.of('binary', 'counts')
        | map(it -> tuple(it, file("${params.base_dir}/raw_matrices/matrix.${it}.mtx.gz")))
        | combine(Channel.fromPath(params.index_file))
        | combine(singletons_strategy)
        | filterAndConvertToNumpy
}

workflow getMasks {
    params.base_dir = params.outdir
    singletons_strategy = Channel.of('keep_all')
    matrices = Channel.fromPath("${params.base_dir}/raw_matrices/matrix.binary.mtx.gz")
        | combine(Channel.fromPath(params.index_file))
        | combine(singletons_strategy)
        | filter_masterlist
}


// DEFUNC

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