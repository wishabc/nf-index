# nf-index
Nextflow pipelines to construct a chromatin accessibility peak index and do follow-up analysis

# Pipelines:
- build_masterlist.nf - Build an index of accessible elements using the approach described in [Meuleman et al](https://www.nature.com/articles/s41586-020-2559-3).
- generate_matrices.nf - Using constructed index as a scaffold to generate count (# of reads overlapping DHS) and binary (absence/presence of a peak) matricies.
- filter_peaks.nf - Filter peaks and convert data to np binary format for follow-up analysis. We filter:<br>
  1) Peaks overlapping ENCODE blacklisted regions
  2) Low-signal singletons (low signal peaks identified in only one sample)
  3) Peaks located on non-autosomal chromosomes (for NMF analysis and normalization)
- normalize_signal.nf - Normalize filtered count matrix by running lowess normalization followed by DEseq2 variance-stabilizing transformation (VST). It can also be used to normalize new samples using existing parameters.
- nmf.nf - DEFUNC, Run NMF using normalized matrix
- variance_partition.nf - Run variance partition using normalized matrix
### Main workflow
main.nf - Annotate index and run `build_masterlist, generate_matrices, filter_peaks and normalize_signal` pipelines.

# Input data
TODO:
