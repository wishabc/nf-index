#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(reticulate)
library(stringr)
library(data.table)
library(DESeq2)
library(SummarizedExperiment)


np <- import("numpy", convert=FALSE)

# Input parameters:
# 1) counts matrix 
# 2) normalization factors - pre-computed normalization factors. 
# If non-existent file, perform usual VST
# 3) sample_names - sample names order, will be used as colnames for the matrices
# 4) metadata_file - metadata, one row for each sample. 
# Should be in the same order as sample_names
# 5) prefix - prefix of output filenames:
# 6) params_file - (optional!) already calculated VST parameters (from previous run
# of this script). Normalize the data according to the provided model.

# The script creates two files <prefix>.params.RDS and <prefix>.dds.RDS

args = commandArgs(trailingOnly=TRUE)

# test if there is at least four arguments: if not, return an error
if (length(args) < 4) {
  stop("At least four input arguments should be supplied", call.=FALSE)
}

sample_names <- fread(args[1], sep="\n", header=FALSE)$V1

sample_meta <- as.data.frame(fread(args[2], stringsAsFactors = TRUE))
row.names(sample_meta) <- sample_meta$ag_id
sample_meta <- sample_meta[sample_names, ]


print("Reading counts")
counts <- np$load(args[3])
print('Converting counts to R')
counts <- as.matrix(py_to_r(counts))
storage.mode(counts) <- "integer"
colnames(counts) <- sample_names

dhs_id <- fread(args[4], sep="\n", header=FALSE)$V4
row.names(counts) <- dhs_id

print(colnames(counts))
print(rownames(sample_meta))
stopifnot(identical(colnames(counts), rownames(sample_meta)))

print('Making DESeq dataset')
dds <- DESeqDataSet(
    SummarizedExperiment(assays=list(counts=counts), colData=sample_meta),
    design=~1
)

rm(counts)
gc()

# Provide NULL or non-existent norm_factors file for conventional VST
norm_factors_path <- NULL 
if (length(args) >= 6) {
  print("Taking pre-computed params")
  norm_factors_path <- args[6]
}


if (is.null(norm_factors_path)) {
  print('Reading norm factors')
  norm_factors <- as.matrix(py_to_r(np$load(norm_factors_path)))
  print('Converting norm factors to R')
  normalizationFactors(dds) <- norm_factors
  rm(norm_factors)
  gc()
} else {
  print("Calculating size factors")
  dds <- estimateSizeFactors(dds)
}

dds_name <- args[5]
timing <- system.time({
  saveRDS(dds, dds_name, compress = FALSE)
})
