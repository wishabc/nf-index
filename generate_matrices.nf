#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { add_matrices_to_anndata; convert_to_numpy; extract_meta_from_anndata } from "./converters"

params.conda = "$moduleDir/environment.yml"


process convert_regions_to_summits {
    conda params.conda

    input:
        path masterlist
    
    output:
        tuple path(summit_masterlist), path(inverse_argsort)
    
    script:
    summit_masterlist = "masterlist.summits.bed"
    inverse_argsort = "masterlist.summits.inverse_argsort.txt"
    """
    awk -F'\t' -v OFS='\t' \
        '{ print \$1, \$5, \$5 + 1, \$4, NR }' ${masterlist} \
        | sort-bed - > ${summit_masterlist}

    cat ${summit_masterlist} \
        | awk '{print NR, \$NF}' \
        | sort -k2,2n | awk '{print \$1 - 1}' > ${inverse_argsort}
    """
}

process restore_masterlist_order {

    conda params.conda
    tag "${prefix}"

    input:
        tuple val(prefix), path(np_array, name: "input/*"), path(inverse_argsort)
    
    output:
        tuple val(prefix), path(name)
    
    script:
    name = np_array.name
    """
    python3 $moduleDir/bin/helpers/restore_masterlist_order.py \
        ${np_array} \
        ${inverse_argsort} \
        ${name}
    """
}


// Processes below create a np array of shape n_dhs x 1
process project_peak_calls {
    conda params.conda
    tag "${sample_id}"
    scratch true

    input:
		tuple path(masterlist), val(sample_id), path(peaks_file)
    
    output:
        tuple val(suffix), path(name)
    
    script:
    suffix = 'binary'
    name = "${sample_id}.${suffix}.npy"
    """
    zcat ${peaks_file} \
        | grep -v "^#" \
        | bedmap --fraction-either 0.5 \
        --indicator \
        --sweep-all \
        ${masterlist} \
        - > tmp.txt

    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype bool
    """
}


process extract_max_density {
    conda params.conda
    tag "${sample_id}"
    scratch true

    input:
        tuple path(masterlist), val(sample_id), path(density_bw)
    
    output:
        tuple val(suffix), path(name)
    
    script:
    suffix = "density"
    name = "${sample_id}.${suffix}.npy"
    """
    bigWigToBedGraph ${density_bw} tmp.bg

    cat tmp.bg \
        | awk -v OFS='\t' '{print \$1,\$2,\$3,"${sample_id}",\$4}' \
        | bedmap --sweep-all \
            --delim "\t" \
            --max ${masterlist} - \
        | sed 's/\\<NAN\\>/0/g' \
        > tmp.txt

    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype float
    """
}

process count_tags {
	tag "${sample_id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple path(saf), val(sample_id), path(bam_file), path(bam_file_index)

	output:
		tuple val(suffix), path(name)

	script:
    suffix = "counts"
    name = "${sample_id}.${suffix}.npy"
    ext = bam_file.extension
	"""
    if [ ${ext} != 'bam' ]; then 
        samtools view -bh \
            --reference ${params.genome_fasta} \
            ${bam_file} > align.bam
        samtools index align.bam
    else
        ln -s ${bam_file} align.bam
        ln -s ${bam_file_index} align.bam.bai
    fi

    count=`samtools stats align.bam \
        | grep "^SN" \
        | grep "mapped and paired" \
        | cut -f3` && [[ \$count -gt 0 ]] && tag="-p" || tag=""

	featureCounts -a ${saf} -O -o counts.txt -F SAF \$tag align.bam
	cat counts.txt \
        | awk 'NR > 2 {print \$(NF)}' > tmp.txt
        
    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype int
	"""
}

process extract_max_pval {
    tag "${sample_id}"
    conda "/home/sabramov/miniconda3/envs/hotspot3"

    input:
        tuple path(masterlist), val(sample_id), path(pvals_parquet)

    output:
        tuple val(suffix), path(name)

    script:
    suffix = 'neglog10_pvals'
    name = "${sample_id}.${suffix}.npy"
    """
    hotspot3-pvals \
        ${pvals_parquet} \
        ${masterlist} \
        ${name} \
        --chrom_sizes ${params.nuclear_chrom_sizes} \
        --format npy
    """
}

// returns a np array of shape n_dhs x 2
process extract_bg_mean {
    conda params.conda
    tag "${sample_id}"
    label "bg_params"

    input:
        tuple path(masterlist), val(sample_id), path(bg_params_tabix)
    
    output:
        tuple val(suffix), path(name)
    
    script:
    suffix = "mean_bg_agg_cutcounts"
    name = "${sample_id}.${suffix}.npy"
    """
    echo -e "chrom_masterlist\tstart_masterlist\tend_masterlist\tdhs_id\tsummit\t\$(head -1 <(zcat ${bg_params_tabix}))" > tmp.bed

    zcat ${bg_params_tabix} \
        | grep -v "^#" \
        | grep segment \
        | bedtools intersect \
            -loj \
            -a ${masterlist} \
            -b stdin \
            -sorted \
            >> tmp.bed

    python3 $moduleDir/bin/helpers/extract_bg_params.py tmp.bed ${name}
    """
}

// merge a np array of shape n_dhs x 2
process generate_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "${name}"
	
	label "highmem"
    conda "${params.conda}"
	scratch true
    tag "${prefix}"

	input:
		tuple val(prefix), path(files), path(samples_order)

	output:
		tuple val(prefix), path(name)

	script:
    name = "matrix.${prefix}.npy"
	"""
    python3 $moduleDir/bin/helpers/matrix_from_vectors.py \
        ${prefix} \
        ${samples_order} \
        ${name} \
        --input_ext npy
	"""
}


workflow generateMatrices {
    take:
        data // masterlist, samples_order, saf_masterlist, id, bam, bam_index, peaks, density, fit_stats, hotspot3_pvals_parquet
    main:

        binary_cols = data
            | map(it -> tuple(it[0], it[3], it[6]))
            | project_peak_calls
        
        count_cols = data
            | map(it -> tuple(it[2], it[3], it[4], it[5]))
			| count_tags
        
        max_pvals = data
            | map(it -> tuple(it[0], it[3], it[9]))
            | extract_max_pval

        // Summit based extraction
        summits_masterlist = data
            | map(it -> it[0])
            | convert_regions_to_summits // summits_masterlist, inverse_argsort

        density_cols = data
            | map(it -> tuple(it[3], it[7])) // samples_order, density_bw
            | combine(summits_masterlist.map(it -> it[0])) // summits_masterlist
            | map(it -> tuple(it[3], it[0], it[1], it[2]))
            | extract_max_density

        summit_based_cols = data
            | map(it -> tuple(it[3], it[8]))
            | combine(summits_masterlist.map(it -> it[0])) // summits_masterlist
            | map(it -> tuple(it[2], it[0], it[1]))
            | extract_bg_mean
            | mix(density_cols) 
            | combine(summits_masterlist.map(it -> it[1])) 
            | restore_masterlist_order
        

        samples_order = data
            | map(it -> it[1])
            | first()

        out = binary_cols
            | mix(count_cols)
            | mix(max_pvals)
            | mix(summit_based_cols)

            | combine(samples_order.countLines().toInteger()) // suffix, column, samples_order
            | map(it -> tuple(groupKey(it[0], it[2]), it[1]))
            | groupTuple(by: 0) // suffix, columns
            | combine(samples_order) // suffix, columns, samples_order
            | generate_matrix
    emit:
        out
}


workflow generateMatricesFromAnndata {
    take:
        bams_hotspots
    main:
        print "Looking for anndata input (generated by build_index.nf) in - params.index_anndata = ${params.index_anndata}"

        index_anndata = Channel.of(params.index_anndata)

        matrices = index_anndata
            | extract_meta_from_anndata // masterlist , saf_masterlist, samples_order, samples_meta
            | map(it -> tuple(it[0], it[2], it[1])) // masterlist, samples_order, saf_masterlist
            | combine(bams_hotspots)
            | generateMatrices
            | map(it -> it[1])
            | collect(sort: true, flat: true)

        add_matrices_to_anndata(
            index_anndata,
            matrices
        )
    emit:
        add_matrices_to_anndata.out
}

// main workflow 
workflow {
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row[params.matrix_peaks_column]),
            file(row.normalized_density_bw),
            file(row.hotspot3_fit_stats_file),
            file(row.hotspot3_pvals_parquet)
        ))
        | generateMatricesFromAnndata
}


workflow extractDensityAtSummit {
    println "Extracting normalized density at summits from ${params.index_anndata}"
    data = Channel.of(params.index_anndata)
        | extract_meta_from_anndata // masterlist , saf_masterlist, samples_order, samples_meta
    
    summits_masterlist = data
        | map(it -> it[0])
        | convert_regions_to_summits // summits_masterlist, inverse_argsort

    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.normalized_density_bw)
        ))
        | combine(summits_masterlist.map(it -> it[0])) // summits_masterlist
        | map(it -> tuple(it[2], it[0], it[1]))
        | extract_max_density
        | combine(summits_masterlist.map(it -> it[1])) // inverse_argsort
        | restore_masterlist_order
        | groupTuple()
        | combine(data.map(it -> it[2])) // suffix, column, samples_order
        | generate_matrix
}

////////////////////// DEFUNC  /////////////
workflow extractBGParams {
    println "Extracting bg params at ${params.reference_bed}"
    reference_bed = Channel.fromPath(params.reference_bed)
    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.peak_stats),
            file("${row.peak_stats}.tbi")
        ))
        | combine(reference_bed)
        | map(it -> tuple(it[3], it[0], it[1], it[2]))
        | extract_bg_params
}
