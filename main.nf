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


workflow generateMatrix {
	take:
		BAMS_HOTSPOTS
	main:
		COUNT_FILES = count_tags(BAMS_HOTSPOTS)

		collected_files = COUNT_FILES.toList().transpose().toList()
		collected_files.view()
		generate_count_matrix(collected_files)
	emit:
		generate_count_matrix.out
}
workflow {
	BAMS_HOTSPOTS = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.ag_id, row.bam_file, row.hotspots_file) }
	generateMatrix(BAMS_HOTSPOTS)
}