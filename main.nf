#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process count_tags {
	tag "${indiv_id}"
	conda params.conda

	input:
		tuple val(indiv_id), path(bam_file), path(bam_file_index), path(peaks_file)

	output:
		tuple val(indiv_id), path("${prefix}.counts.txt"), path("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	bedtools multicov -bams ${bam_file} \
		-bed ${params.index_file} \
		| awk '{print \$(NF)}' > ${prefix}.counts.txt
	bedmap --indicator --sweep-all ${params.index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "matrix.all*"
	publishDir params.outdir, pattern: "indivs_order.txt"

	input:
		tuple val(indiv_ids), path(count_files), path(bin_files)

	output:
		tuple path("matrix.all.signal.txt.gz"), path("matrix.all.peaks.txt.gz"), emit: matrices
		path "indivs_order.txt", emit: indivs

	script:
	indiv_ids_join = indiv_ids.join("\t")
	"""
	echo "${indiv_ids_join}" > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}

process filter_index {
	publishDir "${params.outdir}"
	scratch true
	conda params.conda

	output:
		path filtered_mask, emit: mask
		path filtered_index, emit: filtered_index
		

	script:
	filtered_index = "masterlist.filtered.bed"
	filtered_mask = 'filtered_peaks_mask.txt'
	"""
	bedmap --indicator --sweep-all --bp-ovr 1 ${params.index_file} \
        ${params.encode_blacklist_regions} > blacklisted_mask.txt
	
	python3 $moduleDir/bin/filter_index.py \
        ${params.index_file} \
        blacklisted_mask.txt \
        ${filtered_mask}
		${filtered_index}
	"""
}

process apply_filter_to_matrix {
	publishDir "${params.outdir}"
	label "bigmem"

	input:
		tuple path(signal_matrix), path(peaks_matrix)
		path filtered_mask
	
	output:
		tuple path(signal_filt_matrix), path(peaks_filt_matrix)
	
	script:
	signal_filt_matrix = "signal.filtered.matrix.npy"
	peaks_filt_matrix = "binary.filtered.matrix.npy"
	"""
	(
		trap 'kill 0' SIGINT; \
		
		python3 $moduleDir/bin/convert_to_numpy.py \
			${signal_matrix} \
			${signal_filt_matrix} \
			--mask ${filtered_mask} & \
		
		python3 $moduleDir/bin/convert_to_numpy.py \
			${peaks_matrix} \
			${peaks_filt_matrix} \
			--mask ${filtered_mask} & \
		
		wait \
	)
	"""
}


process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.normed.npy"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.params.npz"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.scale_factors.npy"

	input:
		tuple path(signal_matrix), path(peaks_matrix)

	output:
		tuple path(signal_matrix), path("${prefix}.scale_factors.npy"), emit: scale_factors
		path "${prefix}.normed.npy", emit: normed_matrix
		path "${prefix}.params.npz", emit: model_params
		

	script:
	prefix = 'normalized'
	"""
	python3 $moduleDir/bin/lowess.py \
		${peaks_matrix} \
		${signal_matrix} \
		./ \
		--jobs ${task.cpus} \
		--prefix ${prefix}
	"""
}

process reorder_meta {
	conda params.conda
	publishDir "${params.outdir}"

	input:
		path indivs_order

	output:
		path name
	
	script:
	name = "reordered_meta.txt"
	"""
	tr '\t' '\n' < ${indivs_order} > indivs_order_col.txt
	awk -F'\t' 'FNR == NR { lineno[\$1] = NR; next} {print lineno[\$1], \$0;}' indivs_order_col.txt ${params.samples_file} | sort -k 1,1n | cut -d' ' -f2- > ${name}
	"""
}

process deseq2 {
	conda params.conda
	publishDir "${params.outdir}", pattern: "${prefix}*.npy"
	publishDir "${params.outdir}/vst_params", pattern: "${prefix}*.RDS"
	label "bigmem"

	input:
		tuple path(signal_matrix), path(scale_factors)
		path indivs_order
		path meta_path

	output:
		path "${prefix}*"

	script:
	prefix = "deseq.normalized"
	"""
	Rscript $moduleDir/bin/deseq2.R \
	  ${signal_matrix} \
	  ${scale_factors} \
	  ${indivs_order} \
	  ${meta_path} \
	  ${prefix}
	"""

}

workflow generateMatrix {
	take:
		bams_hotspots
	main:
		count_matrices = bams_hotspots
			| count_tags
			| toList()
			| transpose()
			| toList()
			| generate_count_matrix

		mask = filter_index().mask
		out = apply_filter_to_matrix(count_matrices.matrices, mask)
	emit:
		out
		count_matrices.indivs
}

workflow normalizeMatrix {
	take:
		matrices
		indivs_order
	main:
		sf = normalize_matrix(matrices).scale_factors
		new_meta = reorder_meta(indivs_order)
		out = deseq2(sf, indivs_order, new_meta)
	emit:
		out
}

workflow generateAndNormalize {
	take:
		bams_hotspots
	main:
		matrices = generateMatrix(bams_hotspots)
		out = normalizeMatrix(matrices[0], matrices[1])
	emit:
		out
}

workflow {
	bams_hotspots = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.uniq_id, file(row.bam_file), file("${row.bam_file}.crai"), file(row.hotspots_file)))
	generateAndNormalize(bams_hotspots)
}

// Debug code, defunc
workflow test2 {
	mats = Channel.of(tuple(
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/matrix.all.signal.txt.gz'),
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/matrix.all.peaks.txt.gz')
		)
	)
	indivs_order = Channel.of('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/indivs_order.txt')
	mask = filter_index().mask
	out = apply_filter_to_matrix(mats, mask)
	normalizeMatrix(out, indivs_order)
}

workflow test {
	bams_hotspots = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.uniq_id, file(row.bam_file), file("${row.bam_file}.crai"), file(row.hotspots_file)))

	count_matrices = bams_hotspots
		| count_tags
		| toList()
		| transpose()
		| toList()
		| generate_count_matrix
}
