// Create binary matrix workflows

process get_chunks_order {
    input:
        path masterlist
    
    output:
        path name

    script:
    name = "chunks_order.txt"
    """
    cut -f4 ${masterlist} > ${name}
    """
}
process write_rows {

    tag "${chunk_file.simpleName}"
    input:
        tuple path(chunk_file), path(chunks_order), path(samples_order)

    output:
        path density, emit: density
        path binary, emit: binary

    script:
    density = "${chunk_file.baseName}.density.txt"
    binary = "${chunk_file.baseName}.binary.txt"
    """
    /net/seq/data/projects/SuperIndex/erynes/masterLists/writeRowsPerChunkForMatrices \
        ${samples_order} \
        ${chunks_order} \
        ${chunk_file} \
        ${density} \
        ${binary}
    """
}

process collect_chunks {
    publishDir "${params.outdir}/raw_matrices"
    tag "${prefix}"

    input:
        tuple val(prefix), path("chunks/*")
    
    output:
        tuple val(prefix), path(matrix)

    script:
    matrix = "matrix.${prefix}.mtx.gz"
    """
    ls chunks/ | wc -l \
        | awk -v dir="chunks/" \
         '{for(i=1;i<=\$1;i++){printf("%s/chunk%04d.${prefix}.txt ",dir,i)}printf("\\n");}' \
        | xargs cat | gzip > ${matrix}
    """
}

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
		tuple path(saf), path(masterlist), val(id), path(bam_file), path(bam_file_index), val(has_paired)

	output:
		path name

	script:
	tag = has_paired ? '-p' : ''
    name = "${id}.counts.txt"
	"""
	samtools view -bh ${bam_file} > align.bam
	samtools index align.bam

	featureCounts -a ${saf} -o counts.txt -F SAF ${tag} align.bam
	cat counts.txt | awk 'NR > 2 {print \$(NF)}' > ${name}
	"""
}

process generate_count_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "${name}"
	
	label "medmem"
	cpus 2
	scratch true

	input:
		tuple path(files), path(samples_order)

	output:
		tuple val(prefix), path(name)

	script:
    prefix = "count"
    name = "matrix.${prefix}.txt.gz"
	"""
    awk '{printf "%s ", \$0".${prefix}.txt"}' ${samples_order} \
        | xargs paste \
        | gzip > ${name}
	"""
}

workflow generateMatrices {
    take:
        unfiltered_masterlist
        samples_order
        peaks_files
        bams_hotspots
    main:
        chunks_order = unfiltered_masterlist
            | get_chunks_order
        
        rows = peaks_files
            | combine(chunks_order)
            | combine(samples_order)
            | write_rows
        
        binary_data = rows.binary
            | collect(sort: true)
            | map(it -> tuple("binary", it))
        
        binary_and_density = rows.density
            | collect(sort: true)
            | map(it -> tuple("density", it))
            | mix(binary_data)
            | collect_chunks
        
        out = unfiltered_masterlist
			| bed2saf
			| combine(bams_hotspots)
			| count_tags
			| collect(sort: true, flat: true)
            | combine(samples_order)
            | generate_count_matrix	
            | mix(binary_and_density)
    emit:
        out
}

workflow {
    unfiltered_masterlist = Channel.fromPath("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/output/masterlist_DHSs_0802.filtered.annotated.bed")

    samples_order = Channel.of(file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/${params.outdir}/samples_order.txt"))
    peaks_files = Channel.fromPath("/net/seq/data2/projects/sabramov/SuperIndex/dnase-wouter-style-matrices/peaks_list.txt")
        | splitCsv(header: false)
        | map(it -> it[0])
    
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.filtered_alignments_bam),
            file(row?.bam_index ?: "${row.filtered_alignments_bam}.crai"),
            file(row.hotspot_peaks_point1per),
            row.paired_aligned && (row.paired_aligned != 0)
        ))
    
    generateMatrices(unfiltered_masterlist, samples_order, peaks_files, bams_hotspots)
    // file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/work/tmp/2a/10c57903166faeb052d56a4ace1a68/no_core.paths.txt")
    // file("/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_23-09-05/work/tmp/93/23361e11b94264bc72ccbac2af1af4/no_any.paths.txt")
}