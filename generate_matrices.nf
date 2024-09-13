#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { get_samples_order } from "./build_masterlist"
//include { convert_to_h5 } from "./variance_partition"
include { normalizeMatrix } from "./normalize_signal"
include { filterAndConvertToNumpy; convert_to_numpy } from "./filter_peaks"
include { add_matrices_to_anndata } from "./convert_to_anndata"


params.conda = "$moduleDir/environment.yml"


def copy_file(filepath) {
    if (params.index_dir != params.outdir) {
        f = file(filepath)
        if (!f.exists()) {
            error "File not found: ${filepath}"
        }
        file(filepath).copyTo("${params.outdir}/${filepath.name}")
    }

}

process bed2saf {
	conda params.conda

	input:
		path masterlist

	output:
		path name

	script:
	name = "masterlist.saf"
	"""
	cat ${masterlist} \
        | grep -v "^#" \
		| awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' > ${name}
	"""
}


process generate_binary_counts {

    conda params.conda
    tag "${id}"

    input:
		tuple path(masterlist), val(id), path(peaks_file)
    
    output:
        path name
    
    script:
    name = "${id}.binary.txt"
    """
    bedmap --fraction-either 0.5 \
        --indicator \
        <(grep -v "^#" ${masterlist}) \
        ${peaks_file} > ${name}
    """
}


process extract_max_density {
    conda params.conda
    publishDir "${params.outdir}/density"
    tag "${ag_id}"
    scratch true

    input:
        tuple path(masterlist), val(ag_id), path(density_bw)
    
    output:
        path name
    
    script:
    name = "${ag_id}.density.txt"
    """
    bigWigToBedGraph ${density_bw} tmp.bg 
    cat tmp.bg \
        | awk -v OFS='\t' '{print \$1,\$2,\$3,"${ag_id}",\$4}' \
        | bedmap --sweep-all \
            --delim "\t" \
            --max <(grep -v "^#" ${masterlist}) - \
        | sed 's/\\<NAN\\>/0/g'
        > ${name}
    """
}

process count_tags {
	tag "${id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple path(saf), val(id), path(bam_file), path(bam_file_index), path(peaks_file)

	output:
		path name

	script:
    name = "${id}.counts.txt"
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
        | awk 'NR > 2 {print \$(NF)}' > ${name}
	"""
}

process generate_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "${name}"
	
	label "medmem"
	cpus 2
	scratch true

	input:
		tuple val(prefix), path(files)
        path samples_order

	output:
		tuple val(prefix), path(name)

	script:
    name = "matrix.${prefix}.mtx.gz"
	"""
    awk '{printf "%s ", \$0".${prefix}.txt"}' ${samples_order} > file_list.txt

    touch concatenated_output_final.txt

    # Loop through the rest of the batches
    xargs -a file_list.txt -n 500 | while read -r batch; do
        paste "concatenated_output_final.txt" <(paste \$batch) | sed 's/^\t//' > "tmp.txt"
        mv tmp.txt concatenated_output_final.txt
    done

    gzip -c concatenated_output_final.txt > ${name}
	"""
}


workflow generateMatrices {
    take:
        unfiltered_masterlist
        samples_order
        bams_hotspots
    main:
        density_cols = unfiltered_masterlist 
            | combine(bams_hotspots) // masterlist, id, bam, bam_index, peaks, density
            | map(it -> tuple(it[0], it[1], it[5]))
            | extract_max_density
            | collect(sort: true)
            | map(it -> tuple("density", it))

        cols = unfiltered_masterlist
			| bed2saf // 
			| combine(bams_hotspots)
            | map()
			| count_tags
			| collect(sort: true)
            | map(it -> tuple("counts", it))

        all_cols = unfiltered_masterlist
            | combine(bams_hotspots) // masterlist, id, bam, bam_index, peaks, density
            | map(it -> tuple(it[0], it[1], it[4]))
            | generate_binary_counts
            | collect(sort: true)
            | map(it -> tuple("binary", it))
            | mix(cols)
            | mix(density_cols)

        out = generate_matrix(all_cols, samples_order)
    emit:
        out
}


workflow {
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row.peaks_file),
            file(row.normalized_density_bw)
        ))
    unfiltered_masterlist_path = "${params.index_dir}/masterlist_DHSs_all_chunks.${params.masterlist_id}.annotated.bed"
    samples_order_path = "${params.index_dir}/samples_order.txt"
    index_anndata_path = "${params.index_dir}/index.anndata.h5ad"

    if (params.index_dir != params.outdir) {
        copy_file(unfiltered_masterlist_path)
        copy_file(samples_order_path)
        copy_file(index_anndata_path)
    }

    unfiltered_masterlist = Channel.fromPath(unfiltered_masterlist_path)

    samples_order = Channel.fromPath(samples_order_path)


    // Generate matrices
    matrices = generateMatrices(
        unfiltered_masterlist,
        samples_order,
        bams_hotspots
    ) 
        | combine(unfiltered_masterlist)
        | convert_to_numpy
        | map(it -> it[1])
        | collect(sort: true, flat: true)
    
    add_matrices_to_anndata(matrices, index_anndata_path)
    
}
