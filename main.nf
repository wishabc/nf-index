#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

conda_env = "$moduleDir/environment.yml"


process count_tags {
	tag "${indiv_id}"

	publishDir params.outdir + '/count_files'

	input:
		tuple val(index_file), val(indiv_id), val(bam_file), val(peaks_file)

	output:
		tuple val(index_file), val(indiv_id), file("${prefix}.counts.txt"), file("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	$moduleDir/bin/count_tags.py ${bam_file} < ${index_file} > ${prefix}.counts.txt
	
	bedmap --indicator ${index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {

	publishDir params.outdir + '/index'

	input:
		tuple val(index_file), val(indiv_ids), file(count_files), file(bin_files)

	output:
		tuple file("matrix.all.signal.txt.gz"), file("matrix.all.peaks.txt.gz"), file("indivs_order.txt")

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