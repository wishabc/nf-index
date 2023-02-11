#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(reticulate)
library(readr)
library(stringr)
library(data.table)
library(DESeq2)
library(umap)

np <- import("numpy")

# Input parameters:
# 1) counts matrix 
# 2) normalization factors - pre-computed normalization factors. 
# If non-existent file, perform usual VST
# 3) sample_names - sample names order, will be used as colnames for the matrices
# 4) metadata_file - metadata, one row for each sample. 
# Should be in the same order as sample_names
# 5) params_file - (optional!) already calculated VST parameters (from previous run
# of this script). Normalize the data according to the provided model.

args = commandArgs(trailingOnly=TRUE)

# test if there is at least four arguments: if not, return an error
if (length(args) < 4) {
  stop("At least four input files should be supplied", call.=FALSE)
}
# Data reformatting
counts <- np$load(args[1])
# Provide non-existent norm_factors file to perform conventional VST
norm_factors <- ifelse(file.exists(args[2]), np$load(args[2]), NULL)
sample_names <- fread(args[3], header=FALSE)
sample_names <- data.table(sample_names)

counts <- as.data.frame(counts, stringsAsFactors = F)
colnames(counts) <- sample_names

metadata <- read_delim(args[4], delim = '\t', col_names=T)
rownames(metadata) <- metadata$uniq_id

params_file = ifelse((length(args) < 5), NULL, args[5])
# Applying DESEQ with norm_factors
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
if (is.null(norm_factors)) {
  normalizationFactors(dds) <- norm_factors
}

dds <- estimateSizeFactors(dds)
if (is.null(params_file)) {
  dds <- estimateDispersions(dds)
  df <- dispersionFunction(dds)
  saveRDS(df, file=paste(args[5], ".vst.params.RDS", sep=''))
  rm(df)
} else {
  df <- readRDS(params_file)
  dispersionFunction(dds) <- df
}
vsd <- vst(dds, blind = F)
rm(dds)
gc()
np$save(paste(args[5], ifelse(is.null(norm_factors), ".no_sf", ""), ".vst.npy", sep=''), np$array(assay(vsd)))
