#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { generateMatrices; get_samples_order } from "./generate_matrices"
//include { convert_to_h5 } from "./variance_partition"
include { normalizeMatrix } from "./normalize_signal"
include { filterAndConvertToNumpy; filter_masterlist } from "./filter_peaks"

params.conda = "$moduleDir/environment.yml"


process convert_to_anndata {

    script:
    """
    
    """
}


workflow {
    params.base_dir = params.outdir
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
            row.ag_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row.peaks_file)
        ))
    
    unfiltered_masterlist = Channel.fromPath(params.index_file)

    samples_order = Channel.fromPath("${params.base_dir}/samples_order.txt")

    // Generate matrices
    raw_matrices = generateMatrices(
        unfiltered_masterlist,
        samples_order,
        Channel.empty(),
        bams_hotspots
    )

    filters_and_matrices = filterAndConvertToNumpy(unfiltered_masterlist, raw_matrices)

    filters_and_matrices[1]
        | filter(it -> it[0] == "binary")
		| map(it -> it[1])
		| combine(
	        filters_and_matrices[0].map(it -> tuple(it[1], it[2]))
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
