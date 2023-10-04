#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { buildIndex } from "./build_masterlist"
include { generateMatrices; get_samples_order } from "./generate_matrices"
include { convert_to_h5 } from "./variance_partition"

params.conda = "$moduleDir/environment.yml"
params.sample_weights = ""


def non_required_arg(value, key) {
    return value ? "${key} ${value}": ""
}

process apply_filter_and_convert_to_np {
	publishDir "${publishDirectory}", pattern: "${name}"
	label "highmem"
    tag "${prefix}"
	conda params.conda

	input:
		tuple val(prefix), path(matrix), path(mask)
	
	output:
		tuple val(prefix), path(name)
	
	script:
	name = "${prefix}.filtered.matrix.npy"
    publishDirectory = prefix.contains("only_autosomes") ?  "${params.outdir}" : "${params.outdir}/annotations"
    dtype = prefix.contains('binary') ? 'bool' : (prefix.contains('signal') ? 'int' : 'float')
	"""
    python3 $moduleDir/bin/convert_to_numpy.py \
        ${matrix} \
        ${name} \
        --dtype ${dtype} \
        --mask ${mask}
	"""
}


process normalize_matrix {
	conda params.conda
	label "bigmem"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.normed.npy"
	publishDir "${params.outdir}/norm", pattern: "${prefix}.scale_factors.npy"
	publishDir "${params.outdir}/params", pattern: "${prefix}.lowess_params.*"


	input:
        path peaks_matrix
		path signal_matrix
		val norm_params

	output:
		tuple path(signal_matrix), path("${prefix}.scale_factors.npy"), emit: scale_factors
		path "${prefix}.normed.npy", emit: normed_matrix
		tuple path("${prefix}.lowess_params.npz"), path("${prefix}.lowess_params.json"), emit: model_params
		

	script:
	prefix = 'normalized.only_autosomes.filtered'
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
	prefix = "deseq_normalized.only_autosomes.filtered"
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


process filter_masterlist {
    conda params.conda

    publishDir "${params.outdir}", pattern: "${filtered_masterlist}"
    publishDir "${params.outdir}/masks", pattern: "*.mask.txt"

    input:
        path masterlist
    
    output:
	    tuple path(name), path(mask), path(filtered_masterlist), path(autosomes_mask)

    script:
    prefix = "masterlist"
    name = "${prefix}_DHSs.blacklistfiltered.bed"
    mask = "${prefix}.bad_dhs.mask.txt"
    autosomes_mask = "${prefix}.filtered.autosomes.mask.txt"
    filtered_masterlist = "${prefix}.only_autosomes.filtered.bed"
    """
    bedmap --bases ${masterlist} ${params.encode_blacklist_regions} \
        |  awk -F'\t' '{ if(\$1 > 0) print (NR-1)}' \
        > blacklist_rows.txt

    python3 $moduleDir/bin/DHS_filter.py \
        ${prefix} \
        .5 \
        blacklist_rows.txt \
        ${masterlist} \
        ${name} \
        ${mask}
    
    cat ${masterlist} \
		| awk '{print (\$1 ~ /^chr[0-9]+/) ? 1 : 0}' \
		> autosomes.mask.txt
    
    awk 'NR==FNR {mask1[NR]=\$0; next} \
        {print mask1[FNR] * \$0}' ${mask} autosomes.mask.txt > ${autosomes_mask}

	awk 'NR==FNR {mask[NR]=\$0; next} mask[FNR] == 1' \
		${autosomes_mask} ${masterlist} > ${filtered_masterlist}
    """
}


process annotate_masterlist {
    conda params.conda
    publishDir "${params.outdir}/annotations"
    scratch true
    label "bigmem"

    input: 
        tuple path(binary_matrix), path(filtered_masterlist), path(mask)

    output:
        path name

    script:
    name = "masterlist_DHSs_${params.masterlist_id}.filtered.annotated.bed"
    """

     echo "${binary_matrix}"
     echo "${filtered_masterlist}"
     head -10 ${filtered_masterlist}

     python $moduleDir/bin/spot1Annotations.py \
        ${binary_matrix} \
	${mask} \
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

workflow normalizeMatrix {
	take:
		matrices
		samples_order
		normalization_params
	main:
        binary_matrix = matrices
            | filter(it -> it[0] == "binary.only_autosomes")
            | map(it -> it[1])
        count_matrix = matrices
            | filter(it -> it[0] == "counts.only_autosomes")
            | map(it -> it[1])
		lowess_params = normalization_params
			| filter { it.name =~ /lowess_params/ }
			| ifEmpty(null)
			| collect(sort: true)

		sf = normalize_matrix(binary_matrix, count_matrix, lowess_params).scale_factors
		deseq_params = normalization_params
			| filter { it.name =~ /params\.RDS/ }
			| ifEmpty(null)
		out = deseq2(sf, samples_order, deseq_params).matrix

        //h5_file = convert_to_h5(binary_matrix, out, samples_order)

	emit:
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
				file(row.hotspot_peaks_point1per)
			))
	emit:
		bams_hotspots
}

workflow {
	bams_hotspots = readSamplesFile()
    // Build index
	index_data = bams_hotspots
		| map(it -> it[3])
		| buildIndex
    
    unfiltered_masterlist = index_data[0]

    samples_order = get_samples_order()
    autosomes_mask = unfiltered_masterlist
        | filter_masterlist // returns filtered_dhs, filtered_dhs_mask, filtered_autosomes_masterlist, filtered_autosomes_mask
        | map(it -> it[3]) // mask
    // Generate matrices
    raw_matrices = generateMatrices(unfiltered_masterlist, samples_order, index_data[1], bams_hotspots)

    autosomes_filtered_matrices = raw_matrices
                | map(it -> tuple("${it[0]}.only_autosomes", it[1]))
                | combine(autosomes_mask)

    matrices = raw_matrices 
        | combine(filter_masterlist.out.map(it -> it[1]))
        | mix(autosomes_filtered_matrices)
        | apply_filter_and_convert_to_np

    // Normalize matrices
    out = normalizeMatrix(matrices, samples_order, Channel.empty())

    // Annotate index
    matrices
        | filter(it -> it[0] == "binary")
		| map(it -> it[1])
		| combine(
	        filter_masterlist.out.map(it -> tuple(it[0], it[1]))
        )
		| annotate_masterlist
}


workflow existingMatrices {
    params.base_dir = "$launchDir/${params.outdir}"

    autosomes_mask = Channel.fromPath(params.index_file)
        | filter_masterlist // returns filtered_dhs, filtered_dhs_mask, filtered_autosomes_masterlist, filtered_autosomes_mask
        | map(it -> it[3]) // mask

    matrices = Channel.of('binary', 'counts')
        | map(it -> tuple("${it}.only_autosomes", file("${params.base_dir}/raw_matrices/matrix.${it}.mtx.gz")))
        | combine(autosomes_mask)
        | apply_filter_and_convert_to_np
    
    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")


    out = normalizeMatrix(matrices, samples_order, Channel.empty())
}


// Debug code below, defunc
workflow annotateMasterlist {
    index_and_mask = Channel.fromPath(params.index_file)
        | filter_masterlist // returns filtered_dhs, filtered_dhs_mask, filtered_autosomes_masterlist, filtered_autosomes_mask
        | map(it -> tuple(it[0], it[1])) // index, mask

    Channel.fromPath("$launchDir/${params.outdir}/annotations/binary.filtered.matrix.npy")
        | combine(index_and_mask)
        | annotate_masterlist
}


workflow npyMatrices {
    params.base_dir = "$launchDir/${params.outdir}"
    matrices = Channel.of('binary.only_autosomes', 'counts.only_autosomes')
        | map(it -> tuple(it, file("${params.base_dir}/${it}.filtered.matrix.npy")))
    
        
    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")

    out = normalizeMatrix(matrices, samples_order, Channel.empty())
}


workflow existingModel {
    params.base_dir = "$launchDir/${params.outdir}"
    matrices = Channel.of('binary.only_autosomes', 'counts.only_autosomes')
        | map(it -> tuple(it, file("${params.base_dir}/${it}.filtered.matrix.npy")))
	
	existing_params = Channel.fromPath("${params.base_dir}/params/*")
		| map(it -> file(it))
    
    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")

    out = normalizeMatrix(matrices, samples_order, existing_params)
}


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
