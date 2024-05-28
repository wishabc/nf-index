# nf-index
Nextflow pipelines to construct a chromatin accessibility peak index and do follow-up analysis

# Pipelines:
- build_masterlist.nf - Build an index of accessible elements using the approach described in [Meuleman et al](https://www.nature.com/articles/s41586-020-2559-3).
- generate_matrices.nf - Using constructed index as a scaffold to generate count (# of reads overlapping DHS) and binary (absence/presence of a peak) matricies.
- filter_peaks.nf - Filter peaks and convert data to np binary format for follow-up analysis. We filter:<br>
  1) Peaks overlapping ENCODE blacklisted regions
  2) Low-signal singletons (low signal peaks identified in only one sample)
  3) Peaks located on non-autosomal chromosomes (for NMF analysis and normalization)
- normalize_signal.nf - Normalize filtered count matrix by running lowess normalization followed by DEseq2 variance-stabilizing transformation (VST). There is a workflow to apply normalization with existing parameters to new samples.
- nmf.nf - Run non-negative matrix factorization (NMF) for set of matrices. More details see below
- variance_partition.nf - Run variance partition using normalized matrix
### Main workflow
main.nf - run `build_masterlist, generate_matrices, filter_peaks and normalize_signal` pipelines + annotate resulting index with genomic annotations.

# Input data
<details open><summary>nmf.nf</summary>

<p>

- **samples_file**: Samples metadata in tsv format. File should contain `id` (unique identifier of the sample) and `sample_label` columns. Other columns are permitted and ignored. Used solely for visualizations.
- **nmf_params_list**: A tsv file with information required to run NMF. Should contain all required columns. NA values in optional columns are permitted. Other, non-specified columns are permitted and ignored. See columns description below:
    + (required) `n_components` - number of components for NMF. 
    + (required) `prefix`: prefix for all input files. n_components will be added to prefix.
    + (required) `matrix_path`: path to matrix to run NMF on. Expected shape: `DHSs x samples`
    + (required) `sample_names`: one-column file without header that contains names of the samples. They should match with values in `id` column of samples metadata (`samples_file` option). Should be a subset of samples defined in `samples_file`.<br> File format: <br>
        <table>
        <tr>
            <td>sample1</td>
        </tr>
        <tr>
            <td>sample2</td>
        </tr>
        <tr>
            <td>...</td>
        </tr>
        <tr>
            <td>sampleX</td>
        </tr>
        </table>
    + (required) `dhs_meta`: metadata for DHSs (rows) in tsv format without header. First 4 columns are treated as `chr`, `start`, `end`, `dhs_id`. `dhs_id` - unique identifier of DHS
    + (optional) `samples_weights`: weights for the samples in tsv format. Useful when you have class imbalance, e.g. more samples of some specific cell type/condition.
    
        Expected to be a two column tsv file: <br>
        <table>
            <tr>
                <th>id</th>
                <th>weight</th>
            </tr>
            <tr>
                <td>Sample1</td>
                <td>0.9</td>
            </tr>
            <tr>
                <td>Sample2</td>
                <td>0.3</td>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
            </tr>
            <tr>
                <td>SampleN</td>
                <td>1.0</td>
            </tr>
        </table>

    + (optional) `peaks_weights`: weights for the DHSs in tsv format. Useful when you have various biases at different peaks. `id` corresponds to dhs_id (4th column in `dhs_meta`)
    
        Expected to be a two column tsv file:<br>
            <table>
        <tr>
            <th>id</th>
            <th>weight</th>
        </tr>
        <tr>
            <td>chunk0001</td>
            <td>0.9</td>
        </tr>
        <tr>
            <td>chunk0002</td>
            <td>0.3</td>
        </tr>
        <tr>
            <td>...</td>
            <td>...</td>
        </tr>
        <tr>
            <td>chunk9999</td>
            <td>1.0</td>
        </tr>
        </table>
</p>
</details>
