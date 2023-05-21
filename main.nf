#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./nmf"

params.conda = "$moduleDir/environment.yml"
params.sample_weights = ""

process bed2saf {
	conda params.conda

	output:
		path name

	script:
	name = "masterlist.saf"
	"""
	cat ${params.index_file} | awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' > ${name}
	"""
}


process count_tags {
	tag "${indiv_id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple val(indiv_id), path(bam_file), path(bam_file_index), path(peaks_file), val(has_paired), path(saf)

	output:
		tuple val(indiv_id), path("${prefix}.counts.txt"), path("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	tag = has_paired ? '-p' : ''
	"""
	samtools view -bh ${bam_file} > align.bam
	samtools index align.bam

	featureCounts -a ${saf} -o counts.txt -F SAF ${tag} align.bam
	cat counts.txt | awk 'NR > 2 {print \$(NF)}' > ${prefix}.counts.txt
	
	bedmap --indicator --sweep-all ${params.index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "matrix.all*"
	publishDir params.outdir, pattern: "indivs_order.txt"
	label "medmem"

	input:

		tuple val(indiv_ids), path(count_files), path(bin_files)

	output:
		tuple path("matrix.all.signal.txt.gz"), path("matrix.all.peaks.txt.gz"), emit: matrices
		path "indivs_order.txt", emit: indivs

	script:
	"""
	cat ${count_files} | xargs -I file basename file | cut -d. -f1 > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}

process filter_index {
	publishDir "${params.outdir}", pattern: "${filtered_index}"
	publishDir "${params.outdir}/masks", pattern: "*_mask.txt"
	scratch true
	conda params.conda

	output:
		path filtered_mask, emit: mask
		path filtered_index, emit: filtered_index
		path 'blacklisted_mask.txt', emit: blacklist_mask
		

	script:
	filtered_index = "masterlist.filtered.bed"
	filtered_mask = 'filtered_peaks_mask.txt'
	"""
	bedmap --indicator --sweep-all --bp-ovr 1 ${params.index_file} \
        ${params.encode_blacklist_regions} > blacklisted_mask.txt
	
	python3 $moduleDir/bin/filter_index.py \
        ${params.index_file} \
        blacklisted_mask.txt \
        ${filtered_mask} \
		${filtered_index}
	"""
}

process apply_filter_to_matrix {
	publishDir "${params.outdir}"
	label "bigmem"
	conda params.conda

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
	publishDir "${params.outdir}/norm", pattern: "${prefix}.scale_factors.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}.params.npz"


	input:
		tuple path(signal_matrix), path(peaks_matrix)
		val norm_params

	output:
		tuple path(signal_matrix), path("${prefix}.scale_factors.npy"), emit: scale_factors
		path "${prefix}.normed.npy", emit: normed_matrix
		tuple path("${prefix}.lowess_params.npz"), path("${prefix}.lowess_params.json"), emit: model_params
		

	script:
	prefix = 'normalized'
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

process reorder_meta {
	conda params.conda

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
	publishDir "${params.outdir}/params", pattern: "${prefix}*.RDS"
	label "bigmem"

	input:
		tuple path(signal_matrix), path(scale_factors)
		path indivs_order
		path meta_path
		val norm_params

	output:
		path "${prefix}*"

	script:
	prefix = "deseq.normalized"
	normalization_params = norm_params ?: ""
	"""
	Rscript $moduleDir/bin/deseq2.R \
	  ${signal_matrix} \
	  ${scale_factors} \
	  ${indivs_order} \
	  ${meta_path} \
	  ${prefix} \
	  ${normalization_params}
	"""
}

workflow generateMatrix {
	take:
		bams_hotspots
	main:
		
		count_matrices = bams_hotspots 
			| combine(bed2saf())
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
		normalization_params
	main:
		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ }
			| ifEmpty(null)
			| collect(sort: true)

		sf = normalize_matrix(matrices, lowess_params).scale_factors
		new_meta = reorder_meta(indivs_order)
		deseq_params = normalization_params
			| filter { it.name =~ /params\.RDS/ }
			| ifEmpty(null)
		out = deseq2(sf, indivs_order, new_meta, deseq_params)
	emit:
		out
}

workflow generateAndNormalize {
	take:
		bams_hotspots
		normalization_params
	main:
		matrices = generateMatrix(bams_hotspots)
		out = normalizeMatrix(matrices[0], matrices[1], normalization_params)
	emit:
		out
}


workflow readSamplesFile {
	main:
		bams_hotspots = Channel.fromPath(params.samples_file)
			| splitCsv(header:true, sep:'\t')
			| map(row -> tuple(
				row.id,
				file(row.filtered_alignments_bam),
				file(row?.bam_index ?: "${row.filtered_alignments_bam}.crai"),
				file(row.hotspot_peaks_point1per),
				row.paired_aligned && (row.paired_aligned != 0)
			))
	emit:
		bams_hotspots
}
workflow {
	bams_hotspots = readSamplesFile()
	out = generateAndNormalize(bams_hotspots, Channel.empty())
}


workflow normalizeExistingMatrices {
	mats = Channel.of(tuple(
		file("$launchDir/${params.outdir}/signal.filtered.matrix.npy"),
		file("$launchDir/${params.outdir}/binary.filtered.matrix.npy")
		)
	)
	indivs_order = Channel.of(file("$launchDir/${params.outdir}/indivs_order.txt"))
	normalizeMatrix(mats, indivs_order, Channel.empty())
	// params.normalization_params_dir = "$launchDir/${params.outdir}/params"
	// file(params.normalization_params_dir, checkIfExists: true, type: 'dir')
	// existing_params = Channel.fromPath("${params.normalization_params_dir}/*")
	// 	| map(it -> file(it))
}

workflow existingModel {
	params.normalization_params_dir = "$launchDir/${params.outdir}/params"
	file(params.normalization_params_dir, checkIfExists: true, type: 'dir')
	existing_params = Channel.fromPath("${params.normalization_params_dir}/*")
		| map(it -> file(it))
	bams_hotspots = readSamplesFile()
	out = generateAndNormalize(bams_hotspots, existing_params)
}



// Debug code below, defunc
workflow test3 {
	mats = Channel.of(tuple(
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/signal.filtered.matrix.npy'),
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/binary.filtered.matrix.npy')
		)
	)
	indivs_order = Channel.of('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/indivs_order.txt')
	params.normalization_params_dir = "/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/params"
	existing_params = Channel.fromPath("${params.normalization_params_dir}/*")
		| map(it -> file(it))
	normalizeMatrix(mats, indivs_order, existing_params)
}
workflow test2 {
	mats = Channel.of(tuple(
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/matrix.all.signal.txt.gz'),
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/matrix.all.peaks.txt.gz')
		)
	)
	indivs_order = Channel.of('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/output/index/indivs_order.txt')
	mask = filter_index().mask
	out = apply_filter_to_matrix(mats, mask)
	normalizeMatrix(out, indivs_order, Channel.empty())
}

workflow test {
	mats = Channel.of(tuple(
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/raw_matrices/matrix.all.signal.txt.gz'),
		file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/raw_matrices/matrix.all.peaks.txt.gz')
		)
	)
	
	mask = filter_index().mask
	out = apply_filter_to_matrix(mats, mask)
}
