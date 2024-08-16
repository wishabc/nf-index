include { collate_and_chunk; process_chunk } from "./build_masterlist"
// Create binary matrix workflows

process get_samples_order {

    publishDir params.outdir
    
    output:
        path name

    script:
    name = "samples_order.txt"
    """
    awk -F"\t" -v col="ag_id" \
        'NR==1{for(i=1;i<=NF;i++)if(\$i==col)c=i}NR>1{if(c)print \$c}' \
            ${params.samples_file} > ${name}
    """
}


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
    scratch true

    input:
        tuple path(chunk_file), path(chunks_order), path(samples_order)

    output:
        path binary

    script:
    binary = "${chunk_file.baseName}.binary.txt"
    """
    $moduleDir/bin/index_scripts/writeRowsPerChunkForMatrices \
        ${samples_order} \
        ${chunks_order} \
        ${chunk_file} \
        tmp.txt \
        ${binary}
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
    bedtools intersect \
        -a <(cut -f1-4 ${masterlist}) \
        -b <(unstarch ${peaks_file} | cit -f1-3) \
        -wo \
        -f 0.5 \
        | sort -k5,5 -k6,6n \
        # choose only one peak in masterlist with LARGEST overlap for each peak in peaks_file
        | awk -F'\t' -v OFS='\t' \
            '{
                overlap = \$NF;
                key = \$5":"\$6"-"\$7;
                if (key != prev_key) {
                    if (NR > 1) {
                        print current_line;
                    }
                    max_overlap = -1;
                    prev_key = key;
                }
                
                if (overlap > max_overlap) {
                    max_overlap = overlap;
                    current_line = \$1;
                }
            } END { print current_line }
            ' \
        | awk -F'\t' \
            'NR==FNR \
                { 
                    if (\$1 in ids) {
                        print "Error: Repetitive element detected in input IDs: " $1 > "/dev/stderr";
                        exit 1;
                    }
                    ids[\$1]; 
                    next 
                
                } \
                { print (\$4 in ids ? 1 : 0) }' \
                - ${masterlist} > ${name}
    """
    // bedmap --fraction-map 0.8 \
    //     --indicator \
    //     ${masterlist} \
    //     ${peaks_file} > ${name}
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
    ls chunks/ \
        | wc -l \
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
		tuple path(saf), path(masterlist), val(id), path(bam_file), path(bam_file_index), path(peaks_file)

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
	cat counts.txt | awk 'NR > 2 {print \$(NF)}' > ${name}
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
        peaks_files
        bams_hotspots
    main:
        cols = unfiltered_masterlist
			| bed2saf
			| combine(bams_hotspots)
			| count_tags
			| collect(sort: true, flat: true)
            | map(it -> tuple("counts", it))

        if (params.method == 'chunks') {
            chunks_order = unfiltered_masterlist
                | get_chunks_order
            
            out = peaks_files
                | combine(chunks_order)
                | combine(samples_order)
                | write_rows
                | collect(sort: true)
                | map(it -> tuple("binary", it))
                | collect_chunks
                | mix(generate_count_matrix(cols, samples_order))

        } else {
            all_cols = unfiltered_masterlist
                | combine(bams_hotspots) // masterlist, id, bam, bam_index, peaks, paired_aligned
                | map(it -> tuple(it[0], it[1], it[4]))
                | generate_binary_counts
                | collect(sort: true)
                | map(it -> tuple("binary", it))
                | mix(cols)

            out = generate_matrix(
                all_cols,
                samples_order
            )
        }

    emit:
        out
}

workflow {
    params.method = "indicator"
    unfiltered_masterlist = Channel.fromPath(params.index_file)
    samples_order = get_samples_order()
    
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row.peaks_file)
        ))
    
    if (params.method == 'chunks') {
        bams_hotspots
            | map(it -> it[3])
            | collect(sort: true)
            | collate_and_chunk
            | flatten()
            | process_chunk
        peaks_files = process_chunk.out[1]
    } else {
        peaks_files = Channel.empty()
    }
    
    generateMatrices(unfiltered_masterlist, samples_order, peaks_files, bams_hotspots)
}
