# nf-index
Nextflow pipelines to build an index of accessible elements and do follow-up analysis

# Requirements
- Nextflow (https://www.nextflow.io/)
- conda (https://conda.io/projects/conda/en/latest/index.html)


# Quick start
- Follow first steps according to [General usage](#general-usage) section.
- Create files required by pipeline (e.g. `samples_meta`, `samples_order` e.t.c.)
- Fill in params required by pipeline in [params.config](#params).
- Run pipeline using the corresponding command as described in [Usage](#usage) section. 
For example, to [run NMF](#nmfnf) use `nextflow run nmf.nf -profile Altius -resume`.

# Description of pipelines:
- build_masterlist.nf - Build an index of accessible elements using the approach described in [Meuleman et al](https://www.nature.com/articles/s41586-020-2559-3).
- generate_matrices.nf - Using constructed index as a scaffold to generate count (# of reads overlapping DHS) and binary (absence/presence of a peak) matricies.
- normalize_signal.nf - Normalize filtered count matrix by running lowess normalization followed by DEseq2 variance-stabilizing transformation (VST). There is a workflow to apply normalization with existing parameters to new samples.
- variance_partition.nf - Run variance partition using normalized matrix

# Usage
### General usage
 0) (Optional) Create conda environment from `environment.yml` file with ```mamba env create -n super-index -f environment.yml```. Activate the environment (conda activate super-index)
 1) Modify `nextflow.config` to computing enviroment specifications
 2) Fill in params paths in ```params.config```. You can also specify parameters in command line. Please find detailed explanation of the parameters in the [Config section](#config).
 3) Run the pipeline with `nextflow run <workflow.nf> -profile Altius -resume`

### build_masterlist.nf
The workflow consists of these processes:
- Collates and chunks all peak files into genomic chunks.
- Builds DHSs per chunk and resolves overlaps (produces “all”, “non-overlapping core”, and “non-overlapping any” versions per chunk).
- Merges chunks into genome-wide masterlists.
- Generates a raw sparse matrix (chunk-wise rows written then concatenated), then converts it to NumPy format.
- Annotates the masterlist (GC content, GENCODE, repeats) and exports an AnnData object.

To run the pipeline, use:
```
nextflow run build_masterlist.nf -profile Altius,new_cluster --resume \
  --samples_file <samples_file> \
  --index_peaks_column <column_name>
```
`--index_peaks_column` default value is `peaks_file_0.001fdr`.

See [build_masterlist.nf params](#build_masterlistnf-params) for the required `<samples_file>` format.

The `<samples_file>` must follow the required name and format as in build_masterlist.nf params

### generate_matrices.nf
The workflow consists of these processes:
- Uses the masterlist (index) as a scaffold to project per-sample files onto DHSs.
- Builds DHS × sample matrices (.npy) for peak overlap and signal/count tracks (e.g., counts, density, hotspot neglog10_pvals, density, mean_bg_agg_cutcounts).
- Writes matrices and adds them back into the index AnnData.

To run the pipeline, use:
```
nextflow run generate_matrices.nf -profile Altius,new_cluster --resume \
  --samples_file <samples_file> \
  --matrix_peaks_column <column_name>
```
Note: Run this workflow from the same launch directory as build_masterlist.nf so the generated index AnnData is detected automatically.

`--matrix_peaks_column` is optional if your peak-path column is `peaks_file_0.01fdr` (default).

See [generate_matrices.nf params](#generate_matricesnf-params) for the required `<samples_file>` format.

### TODO:
Add other workflows description here

# Config
There are two config files in the repository.
- ```nextflow.config``` - contains enviornment configuration. Detailed explanation can be found at https://www.nextflow.io/docs/latest/config.html. 
- ```params.config``` - specifies thresholds and paths to input files.

Parameters for each process can be specified either in ```params.config``` file or with a command line. See below detailed description of parameters for each workflow
## Params
### Common params

- **samples_file**: Samples metadata in tsv format. File should contain `id` (unique identifier of the sample) and `sample_label` columns. Other columns are permitted and ignored.

- **outdir** - directory to save results into. Defaults to `output` folder in the launch directory
- **conda** - (optional) path to installed conda (from environment.yml). If not present, nextflow creates environment from environment.yml (was not tested).


### build_masterlist.nf params

- **samples_file**: A TSV file that must contain the following columns:
  - `sample_id`: unique identifier for each sample
  - a peak-file path column (set by `--index_peaks_column`, default: `peaks_file_0.001fdr`)

  Example format:
            <table>
        <tr>
            <th>sample_id</th>
            <th>peaks_file_0.001fdr</th>
        </tr>
        <tr>
            <td>AG81965</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81965/fdr0.001/AG81965.peaks.fdr0.001.bed.gz</td>
        </tr>
        <tr>
            <td>AG81945</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81945/fdr0.001/AG81945.peaks.fdr0.001.bed.gz</td>
        </tr>
        <tr>
            <td>...</td>
            <td>...</td>
        </tr>
        <tr>
            <td>AG81957</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81957/fdr0.001/AG81957.peaks.fdr0.001.bed.gz</td>
        </tr>
        </table>


### generate_matrices.nf params

- **samples_file**: Required same `sample_id` columns as `build_masterlist.nf` with the additional following columns:
  - a peak-file path column (set by `--matrix_peaks_column`, default: `peaks_file_0.01fdr`).
  - a bigWig file path column `normalized_density_bw`.
  - a hotspot3 stats file path column `hotspot3_fit_stats_file`.
  - a hotspot3 p-values parquet file path column `hotspot3_pvals_parquet`.

  Example additional columns format:
            <table>
        <tr>
            <th>sample_id</th>
            <th>peaks_file_0.01fdr</th>
            <th>normalized_density_bw</th>
            <th>hotspot3_fit_stats_file</th>
            <th>hotspot3_pvals_parquet</th>
        </tr>
        <tr>
            <td>AG81965</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81965/fdr0.001/AG81965.peaks.fdr0.001.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81965/AG81965.normalized_density.bw</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81965/AG81965.fit_stats.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81965/AG81965.pvals.parquet</td>
        </tr>
        <tr>
            <td>AG81945</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81945/fdr0.001/AG81945.peaks.fdr0.001.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81945/AG81945.normalized_density.bw</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81945/AG81945.fit_stats.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81945/AG81945.pvals.parquet</td>
        </tr>
        <tr>
            <td>...</td>
            <td>...</td>
        </tr>
        <tr>
            <td>AG81957</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81957/fdr0.001/AG81957.peaks.fdr0.001.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81957/AG81957.normalized_density.bw</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81957/AG81957.fit_stats.bed.gz</td>
            <td>/net/seq/data2/projects/afathul/klf1_dbd/2026_klf1_dbd/peak_calls/output/AG81957/AG81957.pvals.parquet</td>
        </tr>
        </table>
### TODO: add details about other workflows



