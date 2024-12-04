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

workflow convertToNumpy {
    take:
        raw_matrices

    main:
        raw_np = raw_matrices
            | map(it -> tuple(it[0], it[1]))
            | convert_to_numpy
    emit:
        raw_np
}


workflow {
    params.base_dir = params.outdir
    matrices = Channel.fromPath("${params.base_dir}/raw_matrices/binary.index.raw.matrix.npy")
        | map(it -> tuple("binary", it, file("${params.base_dir}/masterlist_DHSs_all_chunks.Altius.annotated.bed"), 'filter_median'))
        | filter_masterlist
}
