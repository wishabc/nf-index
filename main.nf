#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { non_required_arg } from "./nmf"
include { buildIndex } from "./build_masterlist"

params.conda = "$moduleDir/environment.yml"
params.sample_weights = ""

process bed2saf {
	conda params.conda

	input:
		path masterlist

	output:
		tuple path(name), path(masterlist)

	script:
	name = "masterlist.saf"
	"""
	cat ${masterlist} \
		| awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' > ${name}
	"""
}


process count_tags {
	tag "${id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple path(saf), path(masterlist), val(id), path(bam_file), path(bam_file_index), path(peaks_file), val(has_paired)

	output:
		tuple path("${id}.counts.txt"), path("${id}.bin.txt")

	script:
	tag = has_paired ? '-p' : ''
	"""
	samtools view -bh ${bam_file} > align.bam
	samtools index align.bam

	featureCounts -a ${saf} -o counts.txt -F SAF ${tag} align.bam
	cat counts.txt | awk 'NR > 2 {print \$(NF)}' > ${id}.counts.txt
	
	bedmap --indicator --sweep-all \
		${masterlist} ${peaks_file} > ${id}.bin.txt
	"""
}

process generate_count_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "matrix.all*"
	
	label "medmem"
	cpus 2
	scratch true

	input:
		path files
		path samples_order

	output:
		tuple path("matrix.all.signal.txt.gz"), path("matrix.all.peaks.txt.gz")

	script:
	"""
	(
		trap 'kill 0' SIGINT; \
		awk '{printf "%s ", \$0".counts.txt"}' ${samples_order} \
			| xargs paste \
			| gzip -c > matrix.all.signal.txt.gz & \
		awk '{printf "%s ", \$0".bin.txt"}' ${samples_order} \
			| xargs paste \
			| gzip -c > matrix.all.peaks.txt.gz & \
		wait \
	)
	"""
}

process filter_and_convert_to_np {
	publishDir "${params.outdir}"
	label "bigmem"
	conda params.conda

	input:
		tuple path(signal_matrix), path(peaks_matrix), path(masterlist_file)
	
	output:
		tuple path(signal_filt_matrix), path(peaks_filt_matrix), path(name)
	
	script:
	signal_filt_matrix = "signal.filtered.matrix.npy"
	peaks_filt_matrix = "binary.filtered.matrix.npy"
	name = "masterlist.filtered.bed"
	"""
	# create a mask
	cat ${masterlist_file} \
		| awk '{print (\$1 ~ /^chr[0-9]+/) ? 1 : 0}' \
		> mask.txt
	awk 'NR==FNR {mask[NR]=\$0; next} mask[FNR] == 1' \
		mask.txt ${masterlist_file} > ${name}

	(
		trap 'kill 0' SIGINT; \
		
		python3 $moduleDir/bin/convert_to_numpy.py \
			${signal_matrix} \
			${signal_filt_matrix} \
			--dtype int \
			--mask mask.txt \
		
		python3 $moduleDir/bin/convert_to_numpy.py \
			${peaks_matrix} \
			${peaks_filt_matrix} \
			--dtype bool \
			--mask mask.txt \
		
		wait \
	)
	"""
}


process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.normed.npy"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.scale_factors.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}.lowess_params.*"


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
	prefix = "deseq.normalized"
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


process annotate_masterlist {
    conda params.conda
    publishDir params.outdir
    scratch true

    input: 
        tuple path(binary_matrix), path(filtered_masterlist)

    output:
        path name

    script:
    name = "masterlist_DHSs_${params.masterlist_id}.filtered.annotated.bed"
    """

     python $moduleDir/bin/spot1Annotations.py \
        ${binary_matrix} \
        ${params.samples_file}
 
    bash $moduleDir/bin/simpleAnnotations.sh \
        ${filtered_masterlist} \
        ${params.encode3} \
        ${params.gencode} \
        ${params.gwas_catalog} \
	${params.repeats} 
    
    bash $moduleDir/bin/gencodeAnnotations.sh \
	${filtered_masterlist} \
	${params.gencode} \
	${params.chromInfo} 


    echo -e "#chr\tstart\tend\tdhs_id\ttotal_signal\tnum_samples\tnum_peaks\tdhs_width\tdhs_summit\tcore_start\tcore_end\tmean_signal" > masterlist_header.txt
    echo -e "is_encode3\tencode3_ovr-fraction\tdist_tss\tgene_name\tnum_gwasCatalog_variants\trepeat_class\trepeat_family\trepeat_name" > simpleAnnotations_header.txt
    echo -e "gene_body\texon_subgroup\tis_coding" > gencodeAnnotations_header.txt
    echo -e "spot1_std\tspot1_min\tspot1_mean\tspot1_median\tspot1_max\tspot1_Q1\tspot1_Q3" > spot1_header.txt
    echo -e 'n_gc\tpercent_gc\tn_mappable' > gc_header.txt

	
    faidx -i nucleotide -b ${filtered_masterlist} ${params.genome_fasta} \
		| awk -v OFS="\t" \
            'NR>1 { 
                total=\$4+\$5+\$6+\$7+\$8;
                cg=\$6+\$7;
                print \$1, \$2-1, \$3, cg, cg/total; }' \
		| bedmap --delim "\t" --echo \
			--bases-uniq - ${params.mappable_file} \
        | cut -f4- \
        > gc_content.txt

    paste masterlist_header.txt simpleAnnotations_header.txt gencodeAnnotations_header.txt spot1_header.txt gc_header.txt > header.txt
    paste ${filtered_masterlist} \
        is_encode3.txt \
        dist_gene.txt \
        gwas_catalog_count.txt \
	repeats.txt \
	gencode_annotations.txt \
        spot1_metrics.tsv \
        gc_content.txt \
        | cat header.txt - > ${name}
 
    """
}


workflow generateMatrix {
	take:
		bams_hotspots
		index_file
		samples_order
	main:
		columns = index_file
			| bed2saf
			| combine(bams_hotspots)
			| count_tags
			| collect(sort: true, flat: true)	
		out = generate_count_matrix(columns, samples_order)
			| combine(index_file)
			| filter_and_convert_to_np
			| map(it -> tuple(it[0], it[1]))
	emit:
		out
		generate_count_matrix.out
}

workflow normalizeMatrix {
	take:
		matrices
		samples_order
		normalization_params
	main:
		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ }
			| ifEmpty(null)
			| collect(sort: true)

		sf = normalize_matrix(matrices, lowess_params).scale_factors
		deseq_params = normalization_params
			| filter { it.name =~ /params\.RDS/ }
			| ifEmpty(null)
		out = deseq2(sf, samples_order, deseq_params).matrix
	emit:
		out
}

workflow generateAndNormalize {
	take:
		data
		index_file
		normalization_params
	main:
		samples_order = data
			| map(it -> it[0])
			| collectFile(
				name: 'samples_order.txt', 
				newLine: true,
				storeDir: params.outdir
			)
		matrices = generateMatrix(data, index_file, samples_order)
		out = normalizeMatrix(matrices[0], samples_order, normalization_params)
	emit:
		matrices[1]
		out
}


workflow readSamplesFile {
	main:
		bams_hotspots = Channel.fromPath(params.samples_file)
			| splitCsv(header:true, sep:'\t')
			| map(row -> tuple(
				row.ag_id,
				file(row.filtered_alignments_bam),
				file(row?.bam_index ?: "${row.filtered_alignments_bam}.crai"),
				file(row.hotspot_peaks_point1per),
				row.paired_aligned && (row.paired_aligned != 0)
			))
	emit:
		bams_hotspots
}

workflow {
	samples = readSamplesFile()
	masterlist = samples
		| map(it -> it[3])
		| buildIndex
	result = generateAndNormalize(
		samples,
		masterlist,
		Channel.empty()
	)
	result[0]
		| map(it -> it[1])
		| combine(masterlist)
		| view()
		| annotate_masterlist
}

workflow existingMasterlist {
	masterlist = Channel.fromPath(params.index_file)
	out = generateAndNormalize(
		readSamplesFile(),
		masterlist,
		Channel.empty()
	)
}

workflow normalizeExistingMatrices {
	mats = Channel.of(tuple(
		file("$launchDir/${params.outdir}/signal.filtered.matrix.npy"),
		file("$launchDir/${params.outdir}/binary.filtered.matrix.npy")
		)
	)
	samples_order = Channel.of(file("$launchDir/${params.outdir}/samples_order.txt"))
	normalizeMatrix(mats, samples_order, Channel.empty())
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
	mats = Channel.of(
		tuple(
			file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/raw_matrices/matrix.all.signal.txt.gz'),
			file('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0108/low_qual_samples/output/raw_matrices/matrix.all.peaks.txt.gz')
		)
	)
	
	mask = filter_index().mask
	out = apply_filter_to_matrix(mats, mask)
}
