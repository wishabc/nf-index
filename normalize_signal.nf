#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { variancePartition } from "./variance_partition"
include { add_normalized_matrices_to_anndata } from "./converters"

params.conda = "$moduleDir/environment.yml"


process extract_from_anndata {
    conda params.conda
    label "highmem"

    input:
        val anndata
    
    output:
        tuple val(prefix), path("matrix.binary.npz"), path("matrix.${params.normalization_layer}.npy"), path(samples_order), path(masterlist), path(samples_meta)

    script:
    prefix = params.dhs_mask_name
    samples_order = "samples_order.txt"
    masterlist = "masterlist.no_header.bed"
    samples_meta = "samples_meta.txt"
    """
    python3 $moduleDir/bin/convert_to_anndata/extract_from_anndata.py \
        ${anndata} \
        ${masterlist} \
        ${samples_order} \
        ${samples_meta} \
        --extra_layers binary ${params.normalization_layer} \
        --dhs_mask_name ${params.dhs_mask_name}
    """
}


process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}", pattern: "${prefix_dir}*.npy"
    publishDir "${params.outdir}/params", pattern: "${prefix_dir}.lowess_params*"
    publishDir "${params.outdir}/qc", pattern: "${prefix_dir}*.pdf"

	input:
        tuple val(prefix), path(peaks_matrix), path(signal_matrix)
		path norm_params, stageAs: "params/*"

	output:
		tuple val(prefix), path("${prefix_dir}.scale_factors.mean_normalized.npy"), emit: scale_factors
        tuple val(prefix), path("${prefix_dir}.log_differences.npy"), emit: log_diffs
		tuple val(prefix), path("${prefix_dir}.lowess_params.npz"), path("${prefix_dir}.lowess_params.json"), emit: lowess_normalization_params
        tuple val(prefix), path("${prefix_dir}.*.pdf"), emit: normalization_qc
		
	script:
    save_dir = 'normalization'
	prefix_dir = "${save_dir}/${prefix}"
	n = norm_params.name != "params/empty.params" ? file(norm_params[0]).baseName : ""
	normalization_params = n ? "--model_params params/${n}" : ""
	"""
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1

	normalize-matrix ${signal_matrix} \
        ${save_dir} \
        --prefix ${prefix} \
		--jobs ${task.cpus} \
        --normalization_type lowess \
        --fit_vs_sample_log_cpm \
		${normalization_params} \

    qc-normalization ${save_dir} \
        --prefix ${prefix} \
        --n_samples 10 \
        --peak_matrix ${peaks_matrix}
	"""
}


process prepare_deseq_dataset {
	conda params.conda
    publishDir "${params.outdir}/deseq2"
    
	label "bigmem"

	input:
		tuple val(prefix), path(scale_factors), path(signal_matrix), path(samples_order), path(masterlist), path(metadata)

	output:
		tuple val(prefix), path(name)

	script:
    name = "${prefix}.deseq_dataset.RDS"
	"""
	Rscript $moduleDir/bin/r_scripts/prepare_deseq_dataset.R \
    	${samples_order} \
        ${metadata} \
		${signal_matrix} \
        ${masterlist} \
        ${name} \
		${scale_factors}
	"""
}


process deseq2_vst {
	conda params.conda
	publishDir "${params.outdir}/normalization", pattern: "${prefix}.vst.npy"
    publishDir "${params.outdir}/params/normalization", pattern: "${prefix}.dispersion_function.RDS"
    publishDir "${params.outdir}/qc", pattern: "${prefix}.dispersions_plot.pdf"
	
    label "bigmem"

	input:
		tuple val(prefix), path(dataset), val(vst_design_formula)
        path vst_params, stageAs: "params/*"

	output:
		tuple val(prefix), path("${prefix}.vst.npy"), emit: vst
        tuple val(prefix), path("${prefix}.dispersion_function.RDS"), emit: dispersion_function
        tuple val(prefix), path("${prefix}.dispersions_plot.pdf"), emit: dispersions_plot

	script:
    p = vst_params.name != "params/empty.params" ? vst_params : ""
	"""
	Rscript $moduleDir/bin/r_scripts/deseq2_vst.R \
        ${prefix} \
		${dataset} \
        '${vst_design_formula}' \
        ${p}
	"""
}


workflow normalizeMatrix {
	take:
		matrices  // prefix, binary_matrix, count_matrix, samples_order, masterlist, samples_file
		normalization_params
	main:
		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ } // TODO modify regex to use prefix
			| ifEmpty(file("empty.params"))
			| collect(sort: true)
        
        vst_params = normalization_params
			| filter { it.name =~ /dispersion_function\.RDS/ } // TODO modify regex to use prefix
			| ifEmpty(file("empty.params"))

		normalization_data = normalize_matrix(
            matrices.map(it -> tuple(it[0], it[1], it[2])),
            lowess_params
        )

        deseq_data = normalization_data.scale_factors // prefix, scale factors
            | join(matrices) // prefix, scale factors, binary_matrix, count_matrix, samples_order, masterlist, samples_file
            | map(it -> tuple(it[0], it[1], it[3], it[4], it[5], it[6])) // prefix, scale factors, count_matrix, samples_order, masterlist, samples_file
            | prepare_deseq_dataset // prefix, dataset
            | combine(Channel.of(params.vst_design_formula)) // prefix, dataset, formula

        vst_data = deseq2_vst(deseq_data, vst_params) // prefix, vst_matrix
        

        variance_partition_data = vst_data.vst
            | join(matrices) // prefix, vst_matrix, binary_matrix, count_matrix, samples_order, masterlist, samples_file
            | map(it -> tuple(it[0], it[1], it[4], it[5], it[6], params.variance_partition_formula))
            | variancePartition // prefix, vp_annotated_masterlist

		out = normalization_data.scale_factors // prefix, scale factors
            | join(matrices.map(it -> tuple(it[0], it[2]))) // prefix, count_matrix
            | join(vst_data.vst)
            | join(normalization_data.lowess_normalization_params)
            | join(vst_data.dispersion_function)
            | join(variance_partition_data) // prefix, scale_factors, vst_matrix, model_params, dispersion_function, vp_annotated_masterlist
            | map(
                it -> tuple(
                    [it[1], it[2], it[3]], // [vst_matrix, scale_factors]
                    [it[4], it[5], it[6]], // [model_params]
                    [it[7]], // bed files
                    [
                        "deseq_design_formula=${params.vst_design_formula}",
                        "variance_partition_formula=${params.variance_partition_formula}",
                        "normalization_layer=${params.normalization_layer}"
                    ] // extras
                )
            )

	emit:
		out // [matrices], [model_params], [bed files], [extras]
}

process extract_normalization_params_from_template {
    conda params.conda
    label "medmem"

    input:
        val anndata

    output:
        path("params/*")

    script:
    """
    mkdir params
    python3 $moduleDir/bin/convert_to_anndata/extract_normalization_params.py \
        ${anndata} \
        params
    """
}

// Re-use existing model from other run (e.g. different samples)
workflow existingModel {
    params.template_anndata = "/path/to/previous/anndata_with_params"
    if (!file(params.template_anndata).exists()) {
        error "Template anndata: ${params.template_anndata} does not exist!"
    } else {
        print "Using existing model from ${params.template_anndata}"
    }
    print "Using layer=${params.normalization_layer} of ${params.matrices_anndata} for normalization"

    existing_params = Channel.of(params.template_anndata)
        | extract_normalization_params_from_template
        | flatten()

    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)

    out = normalizeMatrix(matrices, existing_params)
        
    add_normalized_matrices_to_anndata(anndata, out)
}


workflow existingModelToScaleFactors {
    params.template_anndata = "/path/to/previous/anndata_with_params"
    if (!file(params.template_anndata).exists()) {
        error "Template anndata: ${params.template_anndata} does not exist!"
    } else {
        print "Using existing model from ${params.template_anndata}"
    }
    print "Using layer=${params.normalization_layer} of ${params.matrices_anndata} for normalization"

    existing_params = Channel.of(params.template_anndata)
        | extract_normalization_params_from_template
        | flatten()

    anndata = Channel.of(params.matrices_anndata)
    matrices = extract_from_anndata(anndata)
        | map(it -> tuple(it[0], it[1], it[2]))

    lowess_params = existing_params
        | filter { it.name =~ /lowess_params/ } // TODO modify regex to use prefix
        | collect(sort: true)
    
    normalization_data = normalize_matrix(
        matrices,
        lowess_params
    )
    out = normalization_data.scale_factors
        | join(matrices.map(it -> tuple(it[0], it[2]))) // prefix, count_matrix
        | join(normalization_data.lowess_normalization_params)
        | map(
            it -> tuple(
                [it[1], it[2]], // [scale_factors]
                [it[3], it[4]], // [model_params]
                [], // bed files
                [
                    "normalization_layer=${params.normalization_layer}",
                    "template_anndata=${params.template_anndata}"
                ] // extras
            )
        )
    add_normalized_matrices_to_anndata(anndata, out)
}


// De-novo normalization
workflow {
    print "Using layer=${params.normalization_layer} of ${params.matrices_anndata} for normalization"
    anndata = Channel.of(params.matrices_anndata)
    normalization_input_data = extract_from_anndata(anndata)

    out = normalizeMatrix(normalization_input_data, Channel.empty())

    add_normalized_matrices_to_anndata(anndata, out)
}








process differential_deseq {

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
    Rscript $moduleDir/bin/differential_deseq.R \
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
        | differential_deseq
        | flatten()
        | collectFile(
            storeDir: "${params.outdir}",
            skip: 1,
            keepHeader: true
        ) {
            it -> [ "${it.simpleName}.tsv", it ]
        }
}
