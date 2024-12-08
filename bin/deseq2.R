#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(reticulate)
library(stringr)
library(data.table)
library(DESeq2)

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
if (length(args) < 4) {
  stop("At least four input arguments should be supplied", call.=FALSE)
}

prefix <- args[4]

params_f <- NULL

if (length(args) >= 5) {
  print("Taking existing params")
  params_f <- args[5]
}

print("Reading input matrix")
counts <- np$load(args[1])
# Provide NULL or non-existent norm_factors file for conventional VST
if (is.null(args[2]) | file.exists(args[2])) {
  print('Reading norm factors')
  norm_factors <- np$load(args[2])
} else {
  norm_factors <- NULL
}
sample_names <- fread(args[3], sep="\n", header=FALSE)

# Ensure that sample_names is a vector, not a data table
sample_names <- sample_names$V1

counts <- as.data.frame(counts, stringsAsFactors = F)
colnames(counts) <- sample_names

# metadata <- read_delim(args[4], delim='\t', col_names=T)
# #rownames(metadata) <- metadata$id
colData <- data.frame(row.names=sample_names, sample_ids=sample_names)

print('Making DESeq dataset')
dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~1)
if (is.null(norm_factors)) {
  print("Calculating size factors")
  dds <- estimateSizeFactors(dds)
  suffix <- ".no_sf.vst"
} else {
  print("Applying DESEQ with provided norm_factors")
  normalizationFactors(dds) <- norm_factors
  suffix <- ".sf.vst"
}

if (is.null(params_f)) {
  print('Calculating and saving VST params')
  # code below was copied from https://github.com/mikelove/DESeq2/blob/master/R/vst.R
  # dispersionFunction is not getting saved otherwise
  baseMean <- rowMeans(counts(dds, normalized=TRUE))
  if (sum(baseMean > 5) < 1000) {
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
  dds.sub <- estimateDispersionsFit(dds.sub, fitType="parametric", quiet=TRUE)

  # assign to the full object
  dispersionFunction(dds) <- dispersionFunction(dds.sub)
} else {
  print('Use existing VST params')
  df <- readRDS(params_f)
  dispersionFunction(dds) <- df
}
params_file_name <- paste(prefix, suffix, ".params.RDS", sep='')
if (file.exists(params_file_name)) {
  print(paste('Parameters were not saved. File ', params_file_name, ' exists.', sep=''))
} else {
  saveRDS(dispersionFunction(dds), file=params_file_name)
}

vsd <- varianceStabilizingTransformation(dds, blind = F)
rm(dds)
gc()
np$save(paste(prefix, suffix, ".npy", sep=''), np$array(assay(vsd)))
