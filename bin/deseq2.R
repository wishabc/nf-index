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
# 5) prefix - prefix of output filenames:
# if normalization factors are provided - creates files <prefix>.sf.vst.params.RDS and <prefix>.sf.vst.npy
# else creates files <prefix>.no_sf.vst.params.RDS and <prefix>.no_sf.vst.npy
# 6) params_file - (optional!) already calculated VST parameters (from previous run
# of this script). Normalize the data according to the provided model.

args = commandArgs(trailingOnly=TRUE)

# test if there is at least four arguments: if not, return an error
if (length(args) < 5) {
  stop("At least five input arguments should be supplied", call.=FALSE)
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

prefix <- args[5]
suffix <- ifelse(is.null(norm_factors), ".no_sf.vst", ".sf.vst")

print('Reading params')
params_f <- NULL

if (length(args) < 6) {
  params_f <- args[6]
}

print("Applying DESEQ with norm_factors")
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
if (is.null(norm_factors)) {
  dds <- estimateSizeFactors(dds)
} else {
  normalizationFactors(dds) <- norm_factors
}


if (is.null(params_f)) {
  print('Calculating and saving VST params')
  # code below was copied from https://github.com/mikelove/DESeq2/blob/master/R/vst.R
  # dispersionFunction is not saved otherwise
  baseMean <- rowMeans(counts(dds, normalized=TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, 
  it is recommended to use varianceStabilizingTransformation directly")
  }

  # subset to a specified number of genes with mean normalized count > 5
  dds.sub <- dds[baseMean > 5,]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from=1, to=length(o), length=1000))]
  dds.sub <- dds.sub[idx,]

  # estimate dispersion trend
  dds.sub <- estimateDispersionsGeneEst(dds.sub, quiet=TRUE)
  dds.sub <- estimateDispersionsFit(dds.sub, fitType=fitType, quiet=TRUE)

  # assign to the full object
  dispersionFunction(dds) <- dispersionFunction(dds.sub)
  saveRDS(dispersionFunction(dds), file=paste(prefix, suffix, ".params.RDS", sep=''))
} else {
  print('Use existing VST params')
  df <- readRDS(params_f)
  dispersionFunction(dds) <- df
}
vsd <- vst(dds, blind = F)
rm(dds)
gc()
np$save(paste(prefix, suffix, ".npy", sep=''), np$array(assay(vsd)))
