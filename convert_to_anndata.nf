

process convert_index_to_anndata {
    conda params.conda
    label "highmem"
    publishDir params.outdir

    input:
        tuple path(binary_matrix), path(samples_order), path(masterlist)
        path masks
    
    output:
        path name

    script:
    name = "index.anndata.h5ad"
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
    publishDir params.outdir

    input:
        path index_anndata
        path matrices
    
    output:
        path name

    script:
    name = "index+matrices.anndata.h5ad"
    """
    python3 $moduleDir/bin/convert_to_anndata/matrices_data_to_anndata.py \
        ${index_anndata} \
        ${params.samples_file} \
        ${name} \
        ${matrices}
    """
}
