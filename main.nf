#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { generateMatrices } from "./generate_matrices"
//include { convert_to_h5 } from "./variance_partition"
include { normalizeMatrix } from "./normalize_signal"
include { filterAndConvertToNumpy; filter_masterlist } from "./filter_peaks"

params.conda = "$moduleDir/environment.yml"


def symlink_file(filepath) {
    if (params.index_dir != params.outdir) {
        f = file(filepath)
        file(filepath).copyTo("${params.outdir}/${filepath.name}")
    }

}

workflow {
    // Workflow to generate binary and count matrices from the samples file and existing masterlist
    // Also it filters and creates necessary masks for normalization step

    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row.peaks_file)
        ))
    unfiltered_masterlist_path = "${params.index_dir}/masterlist_DHSs_all_chunks.${params.masterlist_id}.annotated.bed"
    samples_order_path = "${params.index_dir}/samples_order.txt"

    if (params.index_dir != params.outdir) {
        symlink_file(unfiltered_masterlist_path)
        symlink_file(samples_order_path)
    }

    unfiltered_masterlist = Channel.fromPath(unfiltered_masterlist_path)

    samples_order = Channel.fromPath(samples_order_path)


    // Generate matrices
    filters_and_matrices = generateMatrices(
        unfiltered_masterlist,
        samples_order,
        bams_hotspots
    ) 
        | combine(unfiltered_masterlist)
        | filterAndConvertToNumpy
    

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
