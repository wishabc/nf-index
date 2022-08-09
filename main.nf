#!/usr/bin/env nextflow

params.samples_file = '/home/sabramov/projects/SuperIndex/edc-reprocessing/EDC_matrix_master_list.txt'
params.outdir='output'


process count_tags {
	tag "${indiv_id}"

	input:
	tuple val(index_file), val(indiv_id), val(bam_file), val(peaks_file)

	output:
	tuple val(index_file), val(indiv_id), file("${prefix}.counts.txt"), file("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	$baseDir/bin/count_tags.py ${bam_file} ${index_file} > ${prefix}.counts.txt
	
	bedmap --indicator ${index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {

	publishDir params.outdir + '/index'

	input:
	tuple val(index_file), val(indiv_ids), file(count_files), file(bin_files)

	output:
	file "${prefix}.*.txt.gz"

	script:
	"""
	paste  ${count_files} | gzip -c >  matrix.all.signal.txt.gz
	paste - ${bin_files} | gzip -c >  matrix.all.peaks.txt.gz
	"""
}


workflow generateMatrix {
	take:
		BAMS_HOTSPOTS
	main:
		COUNT_FILES = count_tags(BAMS_HOTSPOTS)
		generate_count_matrix(COUNT_FILES.groupTuple(by: 0))
	emit:
		generate_count_matrix.out
}
workflow {
	BAMS_HOTSPOTS = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.index_file, row.ag_id, row.bam_file, row.hotspots_file) }
	generateMatrix(BAMS_HOTSPOTS)
}