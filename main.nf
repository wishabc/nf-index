#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process count_tags {
	tag "${indiv_id}"
	conda "/home/sabramov/miniconda3/envs/babachi"

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
		tuple val(prefix), path(matrix)

	output:
		tuple val(prefix), path(name)
	
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
	label "bigmem"
	publishDir "${params.outdir}/norm"

	input:
		path signal_matrix
		path peaks_matrix

	output:
		path("${prefix}.normed.npy"), emit: normed_matrix
		path("${prefix}.signal.npy"), emit: signal_numpy

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

process get_scale_factors {
	conda params.conda
	publishDir "${params.outdir}"
	label "bigmem"

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
	cpus 10
	label "bigmem"

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
		bams_hotspots
	main:
		count_files = count_tags(bams_hotspots)

		collected_files = count_files.toList().transpose().toList()
		count_matrices = generate_count_matrix(collected_files)
	emit:
		count_matrices.signal
		count_matrices.peaks
		count_matrices.indivs
}

workflow normalizeMatrix {
	take:
		signal_matrix
		peaks_matrix
		indivs_order
	main:
		inp_matrices = signal_matrix
			| map(it -> tuple('signal', it))
			| concat(peaks_matrix.map(it -> tuple('peaks', it)))
			| filter_autosomes
		signal = inp_matrices
			| filter { it[0] == 'signal'}
			| first()
			| map(it -> it[1])
		peaks = inp_matrices
			| filter { it[0] == 'peaks'}
			| first()
			| map(it -> it[1])
		matrices = normalize_matrix(signal, peaks)
		new_meta = reorder_meta(indivs_order)
		signal_np = matrices.signal_numpy
		sf = get_scale_factors(signal_np, matrices.normed_matrix)
		out = deseq2(signal_np, sf, indivs_order, new_meta)
	emit:
		out
}

workflow generateAndNormalize {
	take:
		bams_hotspots
	main:
		matrices = generateMatrix(bams_hotspots)
		out = normalizeMatrix(matrices[0], matrices[1], matrices[2])
	emit:
		out
}

workflow {
	bams_hotspots = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.uniq_id, row.bam_file, row.hotspots_file))
	generateAndNormalize(bams_hotspots)
}


workflow test {
	signal = Channel.of(file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/matrix.all.signal.txt.autosomes.txt.gz'))
	peaks = Channel.of(file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/matrix.all.peaks.txt.autosomes.txt.gz'))
	indivs_order = file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/indivs_order.txt')
	new_meta = file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/reordered_meta.txt')
	
	matrices = normalize_matrix(signal, peaks)
	signal_np = matrices.signal_numpy
	sf = get_scale_factors(signal_np, matrices.normed_matrix)
	out = deseq2(signal_np, sf, indivs_order, new_meta)
}
