
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
        '${params.index_peaks_column}' \
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
        tuple path(matrices), path(normalization_params), path(bed_files), val(extras)
    
    output:
        path(name, type: 'dir')

    script:
    name = "index+matrices+normalized.anndata.zarr"
    extra_args = extras.map(it -> "'" + it + "'").join(' ')
    """
    python3 $moduleDir/bin/convert_to_anndata/add_normalized_matrices_to_anndata.py \
        ${anndata} \
        ${name} \
        --dhs_mask_name ${params.dhs_mask_name} \
        --uns ${extra_args} \
        --varm ${bed_files} \
        --layers ${matrices} \
        --params ${normalization_params}
    """
}

process extract_meta_from_anndata {
    conda params.conda
    label "medmem"

    input:
        val anndata

    output:
        tuple path(masterlist), path(saf_masterlist), path(samples_order), path(samples_meta)

    script:
    masterlist = "masterlist.no_header.bed"
    samples_order = "samples_order.txt"
    saf_masterlist = "masterlist.no_header.saf"
    samples_meta = "samples_meta.txt"
    """
    python3 $moduleDir/bin/convert_to_anndata/extract_from_anndata.py \
        ${anndata} \
        ${masterlist} \
        ${samples_order} \
        ${samples_meta} \
        --matrix_samples_file ${params.samples_file}

    awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' ${masterlist} > ${saf_masterlist}
    """
}
