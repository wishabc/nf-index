#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { convert_to_numpy } from "./filter_peaks"
include { add_matrices_to_anndata } from "./convert_to_anndata"


params.conda = "$moduleDir/environment.yml"


def copy_file(filepath) {
    if (params.index_dir != params.outdir) {
        f = file(filepath)
        if (!f.exists()) {
            error "File not found: ${filepath}"
        }
        file(params.outdir).mkdir()
        f.copyTo("${params.outdir}/${f.name}")
    }

}


process extract_meta_from_anndata {
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
    publishDir "${params.outdir}/density"
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
        | sed 's/\\<NAN\\>/0/g'
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
    xargs -a file_list.txt -n 1000 | while read -r batch; do
        paste "concatenated_output_final.txt" <(paste \$batch) | sed 's/^\t//' > "tmp.txt"
        mv tmp.txt concatenated_output_final.txt
    done

    gzip -c concatenated_output_final.txt > ${name}
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

        all_cols = data
            | map(it -> tuple(it[0], it[3], it[6]))
            | generate_binary_counts
            | mix(cols)
            | mix(density_cols)
            | groupTuple(size: samples_order.countLines().toInteger())

        out = generate_matrix(all_cols, samples_order)
    emit:
        out
}


workflow {
    println "Directory with output of build_masterlist workflow - params.index_dir = ${params.index_dir}"
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
        matrices,
        index_anndata
    )
    
}
