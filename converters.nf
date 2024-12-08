
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

process convert_index_to_anndata {
    conda params.conda
    label "highmem"
    publishDir "${params.outdir}/zarr"

    input:
        tuple path(binary_matrix), path(samples_order), path(masterlist)
    
    output:
        path(name, type: 'dir')

    script:
    name = "index.anndata.zarr"
    """
    python3 $moduleDir/bin/convert_to_anndata/index_data_to_anndata.py \
        ${masterlist} \
        ${samples_order} \
        ${binary_matrix} \
        ${params.samples_file} \
        ${name}
    """
}


process add_matrices_to_anndata {
    conda params.conda
    label "highmem"
    publishDir "${params.outdir}/zarr"

    input:
        val index_anndata
        path matrices
    
    output:
        path(name, type: 'dir')

    script:
    name = "index+matrices.anndata.zarr"
    """
    python3 $moduleDir/bin/convert_to_anndata/matrices_data_to_anndata.py \
        ${index_anndata} \
        ${params.samples_file} \
        ${name} \
        ${matrices}
    """
}

process add_normalized_matrices_to_anndata {
    conda params.conda
    label "highmem"
    publishDir "${params.outdir}"

    input:
        val anndata
        tuple path(matrices), path(normalization_params), path(masterlist_vp), val(formula)
    
    output:
        path(name, type: 'dir')

    script:
    name = "index+matrices+normalized.anndata.zarr"
    """
    python3 $moduleDir/bin/convert_to_anndata/add_normalized_matrices_to_anndata.py \
        ${anndata} \
        ${name} \
        ${masterlist_vp} \
        '${formula}' \
        ${matrices} \
        ${normalization_params}
    """
}