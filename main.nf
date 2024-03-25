#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { buildIndex } from "./build_masterlist"
include { generateMatrices; get_samples_order } from "./generate_matrices"
//include { convert_to_h5 } from "./variance_partition"
include { normalizeMatrix } from "./normalize_signal"
include { filterAndConvertToNumpy, filter_masterlist } from "./filter_peaks"

params.conda = "$moduleDir/environment.yml"
params.sample_weights = ""


process annotate_masterlist {
    conda params.conda
    publishDir "${params.outdir}/annotations"
    scratch true
    label "bigmem"
    errorStrategy 'ignore'

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


    cat ${params.chrom_sizes} \
        | awk -v OFS='\t' '{ print \$1,0,\$2 }' \
        > chrom_sizes.bed


     python $moduleDir/bin/annotations//spot1Annotations.py \
        ${binary_matrix} \
	    ${mask} \
        ${params.samples_file}
 
    bash $moduleDir/bin/annotations/simpleAnnotations.sh \
        ${filtered_masterlist} \
        ${params.encode3} \
        ${params.gencode} \
        ${params.gwas_catalog} \
	    ${params.repeats} 
    
    bash $moduleDir/bin/annotations/gencodeAnnotations.sh \
        ${filtered_masterlist} \
        ${params.gencode} \
        chrom_sizes.bed


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

workflow annotateMasterlist {
    index_and_mask = Channel.fromPath(params.index_file)
        | filter_masterlist // returns filtered_dhs, filtered_dhs_mask, filtered_autosomes_masterlist, filtered_autosomes_mask
        | map(it -> tuple(it[0], it[1])) // index, mask

    Channel.fromPath("${params.outdir}/annotations/binary.filtered.matrix.npy")
        | combine(index_and_mask)
        | annotate_masterlist
}

workflow {
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.filtered_alignments_bam),
            file(row?.bam_index ?: "${row.filtered_alignments_bam}.crai"),
            file(row.hotspot_peaks_point1per)
        ))
    // Build index
	index_data = bams_hotspots
		| map(it -> it[3])
		| buildIndex
    
    unfiltered_masterlist = index_data[0]

    samples_order = get_samples_order()

    // Generate matrices
    raw_matrices = generateMatrices(unfiltered_masterlist, samples_order, index_data[1], bams_hotspots)

    filters_and_matrices = filterAndConvertToNumpy(unfiltered_masterlist, raw_matrices)

    // Normalize matrices
    out = normalizeMatrix(filters_and_matrices[1], samples_order, Channel.empty())

    // Annotate index
    matrices
        | filter(it -> it[0] == "binary")
		| map(it -> it[1])
		| combine(
	        filters_and_matrices[0].map(it -> tuple(it[0], it[1]))
        )
		| annotate_masterlist
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
