

process convert_index_to_anndata {
    conda params.conda
    label "highmem"
    publishDir "${params.outdir}/zarr"

    input:
        tuple path(binary_matrix), path(samples_order), path(masterlist)
        path masks
    
    output:
        path("${name}/**", hidden: true, includeInputs: true)

    script:
    name = "index.anndata.zarr"
    """
    python3 $moduleDir/bin/convert_to_anndata/index_data_to_anndata.py \
        ${masterlist} \
        ${samples_order} \
        ${binary_matrix} \
        ${params.samples_file} \
        ${name} \
        ${masks}
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
        path("${name}/**", hidden: true, includeInputs: true)

    script:
    name = "index+matrices.anndata.zarr"
    """
    echo 1
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
    publishDir "${params.outdir}/zarr"

    input:
        val anndata
        tuple path(matrices), path(normalization_params), path(masterlist_vp), val(formula)
    
    output:
        path("${name}/**", hidden: true, includeInputs: true)

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