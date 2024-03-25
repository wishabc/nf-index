
def non_required_arg(value, key) {
    return value ? "${key} ${value}": ""
}

process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.normed.npy"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.scale_factors.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}.lowess_params.*"


	input:
        path peaks_matrix
		path signal_matrix
		val norm_params

	output:
		tuple path(signal_matrix), path("${prefix}.scale_factors.npy"), emit: scale_factors
		path "${prefix}.normed.npy", emit: normed_matrix
		tuple path("${prefix}.lowess_params.npz"), path("${prefix}.lowess_params.json"), emit: model_params
		

	script:
	prefix = 'normalized.only_autosomes.filtered'
	n = norm_params.size() == 2 ? file(norm_params[0]) : ""
	normalization_params = n ? "--model_params ${n.parent}/${n.baseName}" : ""
	"""
	python3 $moduleDir/bin/lowess.py \
		${peaks_matrix} \
		${signal_matrix} \
		./ \
		--jobs ${task.cpus} \
		--prefix ${prefix} \
		${non_required_arg(params.sample_weights, '--weights')} \
		${normalization_params}
	"""
}

process deseq2 {
	conda params.conda
	publishDir "${params.outdir}", pattern: "${prefix}*.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}*.RDS"
	label "bigmem"

	input:
		tuple path(signal_matrix), path(scale_factors)
		path samples_order
		val norm_params

	output:
		path "${prefix}*.npy", emit: matrix
		path "${prefix}*.RDS", emit: params

	script:
	prefix = "deseq_normalized.only_autosomes.filtered"
	normalization_params = norm_params ?: ""
	"""
	Rscript $moduleDir/bin/deseq2.R \
		${signal_matrix} \
		${scale_factors} \
		${samples_order} \
		${params.samples_file} \
		${prefix} \
		${normalization_params}
	"""
}


workflow normalizeMatrix {
	take:
		matrices
		samples_order
		normalization_params
	main:
        binary_matrix = matrices
            | filter(it -> it[0] == "binary.only_autosomes")
            | map(it -> it[1])
        count_matrix = matrices
            | filter(it -> it[0] == "counts.only_autosomes")
            | map(it -> it[1])
		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ }
			| ifEmpty(null)
			| collect(sort: true)

		sf = normalize_matrix(binary_matrix, count_matrix, lowess_params).scale_factors
		deseq_params = normalization_params
			| filter { it.name =~ /params\.RDS/ }
			| ifEmpty(null)
		out = deseq2(sf, samples_order, deseq_params).matrix

        //h5_file = convert_to_h5(binary_matrix, out, samples_order)

	emit:
		out
}

workflow normalizeNpyMatrices {
    params.base_dir = params.outdir
    matrices = Channel.of('binary.only_autosomes', 'counts.only_autosomes')
        | map(it -> tuple(it, file("${params.base_dir}/${it}.filtered.matrix.npy")))
    
    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")

    out = normalizeMatrix(matrices, samples_order, Channel.empty())
}

workflow existingModel {
    params.base_dir = params.outdir
    matrices = Channel.of('binary.only_autosomes', 'counts.only_autosomes')
        | map(it -> tuple(it, file("${params.base_dir}/${it}.filtered.matrix.npy")))
	
	existing_params = Channel.fromPath("${params.base_dir}/params/*")
		| map(it -> file(it))
    
    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")

    out = normalizeMatrix(matrices, samples_order, existing_params)
}

workflow {

}