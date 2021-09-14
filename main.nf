#!/usr/bin/env nextflow

params.samples_file = '/net/seq/data/projects/regulotyping-h.CD3+/metadata.txt'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'

//Build cell-type specific indicies as well as a combined index
params.build_ct_index = 1 

params.outdir='output'

nuclear_chroms="$params.genome" + ".nuclear.txt"
chrom_sizes="$params.genome"  + ".chrom_sizes.bed"
mappable="$params.genome" + ".K76.mappable_only.bed"
centers="$params.genome" + ".K76.center_sites.n100.nuclear.starch"

// Read samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.indiv_id, row.cell_type, row.bam_file ) }
	.set{ BAMS_HOTSPOTS }


process call_hotspots {
	tag "${indiv_id}:${cell_type}"

	// only publish varw_peaks and hotspots
	publishDir params.outdir + '/hotspots', mode: 'symlink', pattern: "*.starch" 

	module "bedops/2.4.35-typical:modwt/1.0"

	//scratch true

	input:
	file 'nuclear_chroms.txt' from file("${nuclear_chroms}")
	file 'mappable.bed' from file("${mappable}")
	file 'chrom_sizes.bed' from file("${chrom_sizes}")
	file 'centers.starch' from file("${centers}")

	set val(indiv_id), val(cell_type), val(bam_file) from BAMS_HOTSPOTS

	output:
	set val(indiv_id), val(cell_type), val(bam_file), file("${indiv_id}_${cell_type}.varw_peaks.fdr0.001.starch") into PEAKS
	file("${indiv_id}_${cell_type}.hotspots.fdr*.starch")

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

	hotspot2.sh -F 0.05 -f 0.05 -p varWidth_20_${indiv_id}_${cell_type} \
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
		varWidth_20_${indiv_id}_${cell_type} \
		nuclear.cutcounts.starch \
		nuclear.hotspots.fdr0.001.starch \
		../chrom_sizes.bed \
		nuclear.varw_density.fdr0.001.starch \
		nuclear.varw_peaks.fdr0.001.starch \
		\$(cat nuclear.cleavage.total)

		cp nuclear.varw_peaks.fdr0.001.starch ../${indiv_id}_${cell_type}.varw_peaks.fdr0.001.starch
    
    	cp nuclear.hotspots.fdr0.05.starch ../${indiv_id}_${cell_type}.hotspots.fdr0.05.starch
    	cp nuclear.hotspots.fdr0.001.starch ../${indiv_id}_${cell_type}.hotspots.fdr0.001.starch

	rm -rf \${TMPDIR}
	"""
}

PEAKS.into{PEAK_LIST;PEAK_FILES}

PEAK_FILES
	.map{ it -> tuple(it[1], it[3]) }
	.groupTuple(by: 0)
	.tap{PEAK_FILES_BY_CELLTYPE}
	.map{ it -> tuple("all", it[1].flatten()) }
	.set{PEAK_FILES_ALL}

// Include cell-type specific indices or just one large index covering all samples
PEAK_INDEX_FILES = params.build_ct_index ? PEAK_FILES_ALL.concat(PEAK_FILES_BY_CELLTYPE) : PEAK_FILES_ALL

process build_index {
	tag "${cell_type}"
	
	module "R/4.0.5"
	module "bedops/2.4.35-typical"

	publishDir params.outdir + '/index', mode: 'symlink' 

	input:
	set val(cell_type), file('*') from PEAK_INDEX_FILES
	file chrom_sizes from file("${chrom_sizes}")

	output:
	file "masterlist*"
	set val(cell_type), file("masterlist_DHSs_*_nonovl_core_chunkIDs.bed") into INDEX_FILES, INDEX_FILES_FOR_ANNOTATION

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

PEAK_LIST
	.map{ it -> tuple(it[1], it[0], it[2], it[3])}
	.tap{ PEAK_LIST_BY_CELLTYPE }
	.map{ it -> tuple("all", it[1], it[2], it[3])}
	.set{ PEAK_LIST_ALL }

PEAK_LIST_COMBINED = params.build_ct_index ? PEAK_LIST_ALL.concat(PEAK_LIST_BY_CELLTYPE) : PEAK_LIST_ALL

process count_tags {
	tag "${indiv}:${cell_type}"

	//module "python/3.6.4"
	conda '/home/jvierstra/.local/miniconda3/envs/py3.9_default'

	input:
	set val(cell_type), val(indiv_id), val(bam_file), file(peaks_file), file(index_file) from PEAK_LIST_COMBINED.combine(INDEX_FILES, by: 0)

	output:
	set val(cell_type), val(indiv_id), file("${indiv_id}_${cell_type}.counts.txt"), file("${indiv_id}_${cell_type}.bin.txt") into COUNTS_FILES

	script:
	"""
	count_tags.py ${bam_file} < ${index_file} > ${indiv_id}_${cell_type}.counts.txt
	
	bedmap --indicator ${index_file} ${peaks_file} > ${indiv_id}_${cell_type}.bin.txt

	"""
}

process generate_count_matrix {
	tag "${cell_type}"

	publishDir params.outdir + '/index', mode: 'symlink' 

	input:
	set val(cell_type), val(indiv_ids), file(count_files), file(bin_files), file(index_file) from COUNTS_FILES.groupTuple(by: 0).combine(INDEX_FILES_FOR_ANNOTATION, by: 0)

	output:
	file "matrix_*.txt.gz"

	script:
	"""
	echo -n "region_id" > header.txt
	echo -e "\\t${indiv_ids.join("\t")}" >> header.txt

	cat header.txt <(cut -f4 ${index_file} | paste - ${count_files}) | gzip -c >  matrix_counts.${cell_type}.txt.gz
	cat header.txt <(cut -f4 ${index_file} | paste - ${bin_files}) | gzip -c >  matrix_bin.${cell_type}.txt.gz
	"""
}
