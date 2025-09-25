
process core_set {

    conda params.conda
    label 'medmem'
    publishDir "${params.outdir}/core_sets/${params.grouping_column}.${fdr}"
    tag "${prefix}"


    input:
        tuple val(grouping_key), path(pvals), path(anndata), val(fdr)
    
    output:
        tuple val(grouping_key), val(fdr), path("${prefix}.core_set.bed"), path("${prefix}.core_set.npy"), path("${prefix}.saturation_curve.npy"), path("${prefix}.saturation_curve_core.npy"),  path("${prefix}.step_added.npy"), path("${prefix}.mcv_by_step_stats.npy")
    
    script:
    prefix = "${grouping_key}.fdr${fdr}"
    name = "${prefix}.core_set.bed"
    npy_indicator = "${prefix}.npy"
    """
    python3 $moduleDir/bin/core_sets/core_set.py \
        ${prefix} \
        ${params.samples_file} \
        ${params.grouping_column} \
        '${grouping_key}' \
        ${anndata} \
        ${pvals} \
        ${fdr} \
        --method stouffer
    """
}

// FIXME to work with anndata
workflow {
    core_set_fdrs = Channel.from(0.1, 0.05, 0.01, 0.001, 0.0001)
    println "Using ${params.grouping_column} as grouping column"
    params.pvals_matrix = "${params.outdir}/raw_matrices/matrix.neglog10_pvals.npy"
    data = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row[params.grouping_column]
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