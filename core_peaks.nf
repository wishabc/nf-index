include { extract_meta_from_anndata } from "./converters"

process get_mapping {

    conda params.conda
    publishDir "${params.outdir}/core_sets/"
    
    input:
        tuple path(samples_meta), val(grouping_column)
    
    output:
        tuple path(mapping), val(fdrs)
    
    script:
    mapping = "core_set_mapping.${grouping_column}.tsv"
    """
    python3 $moduleDir/bin/core_sets/get_mapping.py \
        ${samples_meta} \
        '${grouping_column}' \
        ${mapping}
    """
}


process core_set {

    conda params.conda
    label 'medmem'
    publishDir "${params.outdir}/core_sets/${fdr}"
    tag "${prefix}"


    input:
        tuple val(grouping_key), val(path_safe_id), val(column), val(fdr), val(anndata)
    
    output:
        tuple val(fdr), val(grouping_key), val(path_safe_id), val(prefix), path("${prefix}.*")
    
    script:
    prefix = "${path_safe_id}.fdr${fdr}"
    """
    python3 $moduleDir/bin/core_sets/core_set.py \
        ${prefix} \
        ${column} \
        '${grouping_key}' \
        ${anndata} \
        --fdr ${fdr} \
        --method ${params.core_set_correction_method}
    """
}

process craft_meta {
    conda params.conda
    publishDir "${params.outdir}/core_sets/"

    input:
        tuple path(mapping), val(fdrs)

    output:
        path name

    script:
    name = 
    fdrs_str = fdrs.join(' ')
    """
    python3 $moduleDir/bin/core_sets/craft_meta.py \
        ${mapping} \
        --fdrs ${fdrs_str}
    """
}

workflow {
    core_set_fdrs = Channel.from(0.1, 0.05, 0.01, 0.001, 0.0001)
    println "Using ${params.grouping_column} as grouping column"

    anndata = Channel.of(params.matrices_anndata)

    mapping = anndata
        | extract_meta_from_anndata // masterlist , saf_masterlist, samples_order, samples_meta
        | map(it -> tuple(it[3], params.grouping_column)) // samples_meta, grouping_column
        | get_mapping
    
    mapping
        | splitCsv(header:true, sep:'\t')
        | map(it -> tuple(it.value, it.path_safe_id, it.column_name))
        | combine(core_set_fdrs) // grouping_key, sample_id, column_name, core_set_fdr
        | combine(anndata)
        | core_set
    
    mapping
        | combine(core_set_fdrs)
        | groupTuple()
        | craft_meta
    
        
}

// FIXME to work with anndata
workflow {
    core_set_fdrs = Channel.from(0.1, 0.05, 0.01, 0.001, 0.0001)
    println "Using ${params.grouping_column} as grouping column"
    
    data = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row[params.grouping_column] // core_ontology_term
            )
        )
        | unique()
        | map(it -> tuple(it[0], file(params.pvals_matrix), file(params.index_anndata)))
        | combine(core_set_fdrs) // grouping_key, pvals, anndata, core_set_fdr
        | core_set
        | collectFile (
            skip: 1,
            keepHeader: true,
            storeDir: "${params.outdir}/core_sets/",
        ) { it -> 
            def parentDir = "${params.outdir}/core_sets/${params.grouping_column}.${it[1]}"
            [ 
                "${params.grouping_column}.core_sets_meta.tsv", 
                "group_key\tfdr\tcore_set_bed\tcore_set_npy\tcore_set_size\tsaturation_curve\tsaturation_curve_core\tstep_added\tmcv_by_step_stats\n${it[0]}\t${it[1]}\t${parentDir}/${it[2].name}\t${parentDir}/${it[3].name}\t${it[2].countLines() - 1}\t${parentDir}/${it[4].name}\t${parentDir}/${it[5].name}\t${parentDir}/${it[6].name}\t${parentDir}/${it[7].name}\n" 
            ] 
        }
}