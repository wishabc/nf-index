#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { convert_to_numpy } from "./filter_peaks"
include { add_matrices_to_anndata } from "./convert_to_anndata"


params.conda = "$moduleDir/environment.yml"


process extract_meta_from_anndata {
    conda params.conda
    label "medmem"

    input:
        path anndata

    output:
        tuple path(masterlist), path(samples_order), path(saf_masterlist)

    script:
    masterlist = "masterlist.no_header.bed"
    samples_order = "samples_order.txt"
    saf_masterlist = "masterlist.no_header.saf"
    """
    python3 $moduleDir/bin/convert_to_anndata/extract_from_anndata.py \
        ${anndata} \
        ${masterlist} \
        ${samples_order}
    awk -v OFS='\t' '{print \$4,\$1,\$2,\$3,"."}' ${masterlist} > ${saf_masterlist}
    """
}


process generate_binary_counts {

    conda params.conda
    tag "${id}"

    input:
		tuple path(masterlist), val(id), path(peaks_file)
    
    output:
        tuple val(suffix), path(name)
    
    script:
    suffix = 'binary'
    name = "${id}.${suffix}.txt"
    """
    bedmap --fraction-either 0.5 \
        --indicator \
        ${masterlist} \
        ${peaks_file} > ${name}
    """
}


process extract_max_density {
    conda params.conda
    tag "${ag_id}"
    scratch true

    input:
        tuple path(masterlist), val(ag_id), path(density_bw)
    
    output:
        tuple val(suffix), path(name)
    
    script:
    suffix = "density"
    name = "${ag_id}.${suffix}.txt"
    """
    bigWigToBedGraph ${density_bw} tmp.bg 
    cat tmp.bg \
        | awk -v OFS='\t' '{print \$1,\$2,\$3,"${ag_id}",\$4}' \
        | bedmap --sweep-all \
            --delim "\t" \
            --max ${masterlist} - \
        | sed 's/\\<NAN\\>/0/g' \
        > ${name}
    """
}

process count_tags {
	tag "${id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple path(saf), val(id), path(bam_file), path(bam_file_index)

	output:
		tuple val(suffix), path(name)

	script:
    suffix = "counts"
    name = "${id}.${suffix}.txt"
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
    tag "${prefix}"

	input:
		tuple val(prefix), path(files), path(samples_order)

	output:
		tuple val(prefix), path(name)

	script:
    name = "matrix.${prefix}.mtx.gz"
	"""
    awk '{print \$0".${prefix}.txt"}' ${samples_order} > file_list.txt

    split -l 400 -d -a 4 \
        file_list.txt \
        batch_file_list_

    flist=\$(ls batch_file_list_* | sort -t '_' -k 2 -n)
    for batch_file in \$flist; do
        batch_output="batch_\${batch_file##*_}.txt"
        paste \$(cat "\$batch_file" | xargs) > "\$batch_output"
    done

    # Combine all batch outputs, sorted numerically, and compress the final result
    paste \$(ls batch*.txt | sort -t '_' -k 2 -n) | gzip -c > ${name}
	"""
}


workflow generateMatrices {
    take:
        data // masterlist, samples_order, saf_masterlist, id, bam, bam_index, peaks, density

    main:
        density_cols = data // masterlist, samples_order, saf_masterlist 
            | map(it -> tuple(it[0], it[3], it[7]))
            | extract_max_density

        cols = data
            | map(it -> tuple(it[2], it[3], it[4], it[5]))
			| count_tags
        
        samples_order = data
            | map(it -> it[1])
            | first()

        out = data
            | map(it -> tuple(it[0], it[3], it[6]))
            | generate_binary_counts
            | mix(cols)
            | mix(density_cols) // suffix, column
            | combine(samples_order.countLines().toInteger()) // suffix, column, samples_order
            | map(it -> tuple(groupKey(it[0], it[2]), it[1]))
            | groupTuple(by: 0) // suffix, columns
            | combine(samples_order) // suffix, columns, samples_order
            | generate_matrix
    emit:
        out
}


workflow generateDensityMatrix {
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.normalized_density_bw)
        ))
    samples_order = Channel.fromPath(params.index_anndata)
        | extract_meta_from_anndata
        | map(it -> it[1])
        | first()

    out = extract_meta_from_anndata.out // masterlist, samples_order, saf_masterlist 
        | map(it -> it[0])
        | combine(bams_hotspots)
        | extract_max_density
        | combine(samples_order.countLines().toInteger()) // suffix, column, samples_order
        | map(it -> tuple(groupKey(it[0], it[2]), it[1]))
        | groupTuple(by: 0) // suffix, columns
        | combine(samples_order) // suffix, columns, samples_order
        | generate_matrix
}


workflow {
    println "Looking for anndata input (generated by build_index.nf) in - params.index_anndata = ${params.index_anndata}"
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row.peaks_file),
            file(row.normalized_density_bw)
        ))

    index_anndata = Channel.fromPath(params.index_anndata)

    matrices = index_anndata
        | extract_meta_from_anndata
        | combine(bams_hotspots)
        | generateMatrices
        | convert_to_numpy
        | map(it -> it[1])
        | collect(sort: true, flat: true)

    
    add_matrices_to_anndata(
        index_anndata,
        matrices,
    )
    
}
