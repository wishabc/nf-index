# nf-index
Nextflow pipelines to build an index of accessible elements and do follow-up analysis

# Requirements
- Nextflow (https://www.nextflow.io/)
- conda (https://conda.io/projects/conda/en/latest/index.html)


# Description of pipelines:
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

# Usage
### General usage
 0) (Optional) Create conda environment from `environment.yml` file with ```conda env create -n super-index -f environment.yml```
 1) Modify `nextflow.config` to computing enviroment specifications
 2) Fill in params paths in ```params.config```. You can also specify parameters in command line. Please find detailed explanation of the parameters in the [Config section](#config).
 3) Run the pipeline with `nextflow run <workflow.nf> -profile Altius -resume`

### NMF.nf
The pipeline consists of two parts:
- Perfroming NMF
- Running QC visualizations

To run both stages of the pipeline use:
```
nextflow run nmf.nf -profile Altius -resume
```

To run just the last, vizualization step (expected to run previous command first):
```
nextflow run nmf.nf -profile Altius -entry visualize --nmf_results_path <launchDir>/output/nmf>
```
The `--nmf_results_path` param can be omitted if you are running the pipeline in the same folder as `nextflow run nmf.nf -profile Altius`. The output files are named according to provided `prefix`. No warning are made in case of name collisions.
### TODO:
Add other workflows description here

# Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```params.config``` - specifies thresholds and paths to input files.

Parameters for each process can be specified either in ```params.config``` file or with a command line. See below detailed description of parameters for each workflow
## Params
<details open><summary>Common parameters</summary>
<p>

- **samples_file**: Samples metadata in tsv format. File should contain `id` (unique identifier of the sample) and `sample_label` columns. Other columns are permitted and ignored.

- **outdir** - directory to save results into. Defaults to `output` folder in the launch directory
- **conda** - (optional) path to installed conda (from environment.yml). If not present, nextflow creates environment from environment.yml (was not tested).

</p>
</details>

<details ><summary>nmf.nf</summary>

<p>

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
    + (required) `dhs_meta`: metadata for DHSs (rows) in tsv format without header. First 4 columns are treated as `chr`, `start`, `end`, `dhs_id`, where `dhs_id` is a unique identifier of a DHS. Other columns are ignored.
    + (optional) `samples_weights`: sample weights in tsv format. NMF prioritizes reconstruction of samples with larger weights. Useful when you have class imbalance, e.g. abundance of samples of some specific cell type/condition.
    
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

    + (optional) `peaks_weights`: weights for the DHSs in tsv format. NMF prioritizes reconstruction of peaks with larger weights. Useful when you have different confidence in different DHSs (rows of the matrix). `id` corresponds to dhs_id (4th column in `dhs_meta`)
    
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

- **dhs_annotations**: (optional, used only for visualizations) A tsv file with DHSs annotations. Should contain `dhs_id` and `dist_tss` columns. Other columns are permitted and ignored. If provided, pipeline plots cumulative distance to tss for DHSs of each component. 

</p>
</details>
TODO: add details about other workflows



