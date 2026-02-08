



process differential_deseq_cell_type {

    memory { 30.GB * task.attempt * task.attempt }
    conda "/home/sabramov/miniconda3/envs/r-jupyter/"
    label "medmem"

    input:
        val chunk

    output:
        path "*${suffix}.tsv"

    script:
    suffix = "deseq_res.${chunk}"
    """
    Rscript $moduleDir/bin/r_scripts/differential_deseq.R \
        ${suffix} \
        ${params.dds} \
        ${params.total_chunks} \
        ${chunk}
    """
}

// nextflow run normalize_signal.nf -profile Altius,new_cluster -entry diffDeseq --dds <path> -resume --total_chunks 1000
workflow diffDeseq {
    params.dds = "/net/seq/data2/projects/sabramov/SuperIndex/hotspot3/w_babachi_new.v23/index/filled_1pr/output/normalization/lowess/ready_for_deseq.reciprocal.mouse+human.normalized.only_autosomes.filtered.no_q.dds.RDS"
    params.total_chunks = 500
    Channel.of(1..params.total_chunks)
        | differential_deseq_cell_type
        | flatten()
        | collectFile(
            storeDir: "${params.outdir}",
            skip: 1,
            keepHeader: true
        ) {
            it -> [ "${it.simpleName}.tsv", it ]
        }
}



process differential_deseq {

    memory { 30.GB * task.attempt * task.attempt }
    conda "/home/sabramov/miniconda3/envs/r-jupyter/"
    label "medmem"
    tag "${meta.prefix}"

    input:
        val meta

    output:
        tuple val(meta), path(name)

    script:
    prefix = meta.prefix
    name = "${prefix}.deseq_results.tsv"
    """
    Rscript $moduleDir/bin/r_scripts/differential_deseq.R \
        ${meta.dataset_anndata} \
        ${prefix} \
        ${name}
    """
}

workflow diffDeseqSpecies {
    Channel.fromPath(params.deseq_meta)
        | splitCsv(header:true, sep:'\t')
        | differential_deseq
        | map(it -> it[1])
        | collectFile(
            storeDir: "${params.outdir}",
            skip: 1,
            keepHeader: true,
            name: "deseq_results.all_cell_types.tsv"
        )
}
