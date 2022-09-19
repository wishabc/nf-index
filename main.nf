#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process count_tags {
	tag "${indiv_id}"
	conda params.conda

	input:
		tuple val(indiv_id), val(bam_file), val(peaks_file)

	output:
		tuple val(indiv_id), file("${prefix}.counts.txt"), file("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	bedtools multicov -bams ${bam_file} -bed ${params.index_file} | awk '{print \$(NF)}' > ${prefix}.counts.txt
	bedmap --indicator ${params.index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {

	publishDir params.outdir + '/index'
	conda params.conda

	input:
		tuple val(indiv_ids), path(count_files), path(bin_files)

	output:
		tuple path("matrix.all.signal.txt.gz"), path("matrix.all.peaks.txt.gz"), path("indivs_order.txt")

	script:
	indiv_ids_join = indiv_ids.join("\t")
	"""
	echo "${indiv_ids_join}" > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}

process normalize_matrix {
	conda params.conda

	cpus params.cpus
	memory params.memory
	publishDir "${params.outdir}/norm"

	input:
		tuple path(signal_matrix), path(peaks_matrix), path(indivs_order)

	output:
		path("${prefix}.normed.npy"), emit: matrix
		path("${prefix}*"), emit: norm_matrices

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
		val metadata_path
		path indivs_order

	output:
		path name
	
	script:
	name = "reordered_meta.txt"
	"""
	while IFS="\t" read LINE
	do
    	grep "\$LINE" "${metadata_path}" >> "${name}"
	done < ${indivs_order} 
	"""
}

process get_scale_factors {
	conda params.conda
	publishDir "${params.outdir}"
	memory params.memory

	input:
		path(signal_matrix)
		path(normed_matrix)

	output:
		path name

	script:
	name = "scale_factors.npy"
	"""
	python3 $moduleDir/bin/get_scale_factors.py \
	  ${signal_matrix} \
	  ${normed_matrix} \
	  ${name}
	"""
}

process deseq2 {
	conda params.conda
	publishDir "${params.outdir}"
	memory params.memory

	input:
		path signal_matrix
		path scale_factors
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
		BAMS_HOTSPOTS
	main:
		COUNT_FILES = count_tags(BAMS_HOTSPOTS)

		collected_files = COUNT_FILES.toList().transpose().toList()
		count_matrices = generate_count_matrix(collected_files)
		
		norm_matrix = normalize_matrix(count_matrices).matrix

		signal_matrix = count_matrices.map(it -> it[0])
		indivs_order = count_matrices.map(it -> it[2])
		new_meta = reorder_meta(params.metadata, indivs_order)

		sf = get_scale_factors(signal_matrix, norm_matrix)

		deseq2(signal_matrix, sf, indivs_order, new_meta)
	emit:
		deseq2.out
}


workflow {
	BAMS_HOTSPOTS = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.ag_id, row.bam_file, row.hotspots_file) }
	generateMatrix(BAMS_HOTSPOTS)
}