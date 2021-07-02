#!/usr/bin/env nextflow

params.outdir='output'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'
params.build_ct_index = 1

nuclear_chroms="$params.genome" + ".nuclear.txt"
chrom_sizes="$params.genome"  + ".chrom_sizes.bed"
mappable="$params.genome" + ".K76.mappable_only.bed"
centers="$params.genome" + ".K76.center_sites.n100.nuclear.starch"

// Read samples file
Channel
	.fromPath("/net/seq/data/projects/regulotyping-h.CD3+/metadata.txt")
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.donor_id, row.cell_type, row.ag_number, row.bamfile ) }
	.into{ SAMPLES_AGGREGATIONS }

SAMPLES_AGGREGATIONS
	.groupTuple(by: [0, 1])
	.map{ it -> tuple(it[0], it[1], it[3].join(" ")) }
	.set{ SAMPLES_AGGREGATIONS_MERGE }

process merge_bamfiles {
	tag "${donor_id}:${cell_type}"

	publishDir params.outdir + '/merged', mode: 'copy' 

	cpus 2

	input:
	set val(donor_id), val(cell_type), val(bam_files) from SAMPLES_AGGREGATIONS_MERGE

	output:
	set val(donor_id), val(cell_type), file('*.bam'), file('*.bam.bai') into BAMS_MERGED 

	script:
	"""
	samtools merge -f -@${task.cpus} ${donor_id}_${cell_type}.bam ${bam_files}
	samtools index ${donor_id}_${cell_type}.bam
	"""
}	

BAMS_MERGED
	.into{ BAMS_MERGED_HOTSPOTS; BAMS_MERGED_COUNTS }

process call_hotspots {
	tag "${donor_id}:${cell_type}"

	publishDir params.outdir + '/hotspots', mode: 'copy' 

	module "bedops/2.4.35-typical:modwt/1.0"

	//scratch true

	input:
	file 'nuclear_chroms.txt' from file("${nuclear_chroms}")
	file 'mappable.bed' from file("${mappable}")
	file 'chrom_sizes.bed' from file("${chrom_sizes}")
	file 'centers.starch' from file("${centers}")

	set val(donor_id), val(cell_type), file(bam_file), file(bam_index_file) from BAMS_MERGED_HOTSPOTS

	output:
	set val(donor_id), val(cell_type), file("${donor_id}_${cell_type}.varw_peaks.fdr0.001.starch") into PEAKS

	script:
	"""

	TMPDIR=\$(mktemp -d)

	samtools view -H ${bam_file} > \${TMPDIR}/header.txt

	cat nuclear_chroms.txt \
	| xargs samtools view -b ${bam_file} \
	| samtools reheader \${TMPDIR}/header.txt - \
	> \${TMPDIR}/nuclear.bam

	PATH=/home/jvierstra/.local/src/hotspot2/bin:\$PATH
	PATH=/home/jvierstra/.local/src/hotspot2/scripts:\$PATH

	hotspot2.sh -F 0.5 -p varWidth_20_default \
		-M mappable.bed \
    	-c chrom_sizes.bed \
    	-C centers.starch \
    	\${TMPDIR}/nuclear.bam \
    	peaks

	cd peaks

	hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch

	rm -f nuclear.varw_peaks.*

	density-peaks.bash \
		\${TMPDIR}\
		varWidth_20_default \
		nuclear.cutcounts.starch \
		nuclear.hotspots.fdr0.001.starch \
		../chrom_sizes.bed \
		nuclear.varw_density.fdr0.001.starch \
		nuclear.varw_peaks.fdr0.001.starch \
		\$(cat nuclear.cleavage.total)

	cp nuclear.varw_peaks.fdr0.001.starch ../${donor_id}_${cell_type}.varw_peaks.fdr0.001.starch

	rm -rf \${TMPDIR}
	"""
}

PEAKS.into{PEAK_LIST;PEAK_FILES}

PEAK_FILES
	.map{ it -> tuple(it[1], it[2]) }
	.groupTuple(by: 0)
	.tap{PEAK_FILES_BY_CELLTYPE}
	.map{ it -> tuple("all", it[1].flatten()) }
	.set{PEAK_FILES_ALL}

// Include cell-type specific indices or just one large index covering all samples
PEAK_INDEX_FILES = params.build_ct_index ? PEAK_FILES_ALL.concat(PEAK_FILES_BY_CELLTYPE) : PEAK_FILES_ALL

process build_index {
	tag "${cell_type}"
	
	module "R/4.0.5"

	publishDir params.outdir + '/index', mode: 'copy' 

	input:
	set val(cell_type), file('*') from PEAK_INDEX_FILES
	file chrom_sizes from file("${chrom_sizes}")

	output:
	set val(cell_type), file ("masterlist*") into INDEX_FILES

	script:
	"""
	ls *.varw_peaks.fdr0.001.starch > filelist.txt

	/home/jvierstra/.local/src/Index/run_sequential.sh \
		\${PWD} \
		${chrom_sizes} \
		filelist.txt \
		${cell_type}

	"""
}

// process count_tags {

// 	input:
// 	set val(donor_id), val(cell_type), file(bam_file), file(bam_index_file) BAMS_MERGED_COUNTS
// 	file 'regions.bed' from INDEX_FILE

// 	output:
// 	set val(donor_id), val(cell_type), file('') into COUNTS_FILES

// 	script:
// 	"""
// 	python3 /home/jvierstra/proj/t_cell_fxn_genotyping/dnase/count_tags_in_region.py \
// 		${bam_file} < regions.bed > ${output_dir}/\${columns[1]}.txt

// 	"""
// }

// process merge_tags {

// }

// process normalize_matrix {


// }
