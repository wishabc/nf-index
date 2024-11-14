nextflow.enable.dsl = 2


params.conda = "$moduleDir/environment.yml"


process convert_to_numpy {
    conda params.conda
    publishDir "${params.outdir}/raw_matrices"
    label "highmem"
    tag "${prefix}"

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


workflow  {
    params.base_dir = params.outdir
    matrices = Channel.fromPath("${params.base_dir}/raw_matrices/binary.index.raw.matrix.npy")
        | map(it -> tuple("binary", it, file("${params.base_dir}/masterlist_DHSs_all_chunks.Altius.annotated.bed"), 'filter_median'))
        | filter_masterlist
}
