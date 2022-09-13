#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process count_tags {
	tag "${indiv_id}"
	conda "${params.conda}"

	input:
		tuple val(indiv_id), val(bam_file), val(peaks_file), path(index_chunk)

	output:
		tuple val(indiv_id), path(index_chunk), path("${prefix}.counts.txt"), path("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	$moduleDir/bin/count_tags.py ${bam_file} < ${index_chunk} > ${prefix}.counts.txt
	
	bedmap --indicator ${index_chunk} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {

	publishDir params.outdir + '/index'
	conda params.conda

	input:
		path(count_files)
		path(bin_files)

	output:
		tuple path("matrix.all.signal.txt.gz"), path("matrix.all.peaks.txt.gz"), path("indivs_order.txt")

	script:
	"""
	cat ${count_files} | xargs -I file basename file | cut -d. -f1 > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}

process split_file {
	scratch true

	input:
		path(input_file)

	output:
		path("out/${prefix}*")

	script:
	prefix = "index_chunk."
	"""
	mkdir out
	cd out
	split -l ${params.chunksize} ${prefix} ${input_file}
	"""
}


workflow generateMatrix {
	take:
		bams_hotspots
		index_chunks
	main:
		count_files = count_tags(bams_hotspots.combine(index_chunks))
		count_indivs = count_files.collectFile(sort: { a,b -> a[1].name <=> b[1].name ?: a[0] <=> b[0] }) {
				item -> ["${item[0]}.count.txt", item[2] + '\n']
		}
		binary_indivs = count_files.collectFile(sort: { a,b -> a[1].name <=> b[1].name ?: a[0] <=> b[0] }) {
				item -> ["${item[0]}.binary.txt", item[3] + '\n']
		}
		generate_count_matrix(count_indivs, binary_indivs)
	emit:
		generate_count_matrix.out
}
workflow {
	bams_hotspots = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple(row.ag_id, row.bam_file, row.hotspots_file))
	index_file = Channel.fromPath(params.index_file)
	generateMatrix(bams_hotspots, split_file(index_file))
}