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

	input:
		tuple val(indiv_ids), path(count_files), path(bin_files)

	output:
		path "matrix.all.signal.txt.gz", emit: signal
		path "matrix.all.peaks.txt.gz", emit: peaks
		path "indivs_order.txt", emit: indivs

	script:
	indiv_ids_join = indiv_ids.join("\t")
	"""
	echo "${indiv_ids_join}" > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}

process filter_autosomes {

	publishDir params.outdir

	input:
		path matrix

	output:
		path name
	
	script:
	name = "${matrix.baseName}.autosomes.txt.gz"
	"""
	len=\$({ cat ${params.index_file} | sed -n '/^chrX/{=;q;}' || true; })
	len=\$((\$len - 1))
	{ zcat ${matrix} | head -n \$len | gzip -c || true; } > ${name}
	"""

}

process normalize_matrix {
	conda params.conda

	cpus params.cpus
	memory { params.memory * task.attempt }
	publishDir "${params.outdir}/norm"

	input:
		path signal_matrix
		path peaks_matrix

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
	tr '\t' '\n' < ${indivs_order} > indivs_order_col.txt
	awk 'FNR == NR { lineno[\$1] = NR; next} {print lineno[\$1], \$0;}' indivs_order_col.txt ${metadata_path} | sort -k 1,1n | cut -d' ' -f2- > ${name}
	"""
}

process get_scale_factors {
	conda params.conda
	publishDir "${params.outdir}"
	memory { params.memory * task.attempt }

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
	memory { params.memory * task.attempt }

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

workflow generateAndNormalize {
	take:
		bams_hotspots
	main:
		generateMatrix(bams_hotspots)
		normalizeMatrix(signal_matrix, peaks_matrix, indivs_order)
}

workflow generateMatrix {
	take:
		bams_hotspots
	main:
		count_files = count_tags(bams_hotspots)

		collected_files = count_files.toList().transpose().toList()
		count_matrices = generate_count_matrix(collected_files)
		signal_matrix = count_matrices.signal
		peaks_matrix = count_matrices.peaks
		indivs_order = count_matrices.indivs
	emit:
		signal_matrix
		peaks_matrix
		indivs_order
}

workflow normalizeMatrix {
	take:
		signal_matrix
		peaks_matrix
		indivs_order
	main:
		matrices = filter_autosomes(Channel.from(signal_matrix, peaks_matrix))
		signal = matrices.first()
		peaks = matrices.last()
		norm_matrix = normalize_matrix(signal, peaks).matrix
		new_meta = reorder_meta(params.metadata, indivs_order)
		sf = get_scale_factors(signal, norm_matrix)
		deseq2(signal, sf, indivs_order, new_meta)
	emit:
		deseq2.out
}

workflow test {
	signal_matrix = file('/net/seq/data/projects/sabramov/SuperIndex/raj+atac_2022-09-10/output/matrix.all.signal.txt.autosomes.txt.gz')

	indivs_order = file('/net/seq/data/projects/sabramov/SuperIndex/raj+atac_2022-09-10/output/indivs_order.txt')
	norm_matrix = file('/net/seq/data/projects/sabramov/SuperIndex/raj+atac_2022-09-10/output/output/matrix_sorted.normed.npy')

	new_meta = file('/net/seq/data/projects/sabramov/SuperIndex/raj+atac_2022-09-10/output/reordered_meta.txt')

	sf = get_scale_factors(signal_matrix, norm_matrix)
	deseq2(signal_matrix, sf, indivs_order, new_meta)
}

workflow {
	BAMS_HOTSPOTS = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.ag_id, row.bam_file, row.hotspots_file) }
	generateMatrix(BAMS_HOTSPOTS)
}