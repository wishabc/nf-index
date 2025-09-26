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

prefix <- args[4]

params_f <- NULL

if (length(args) >= 5) {
  print("Taking existing params")
  params_f <- args[5]
}

print("Reading counts")
counts <- np$load(args[1])

print('Converting counts to R')
counts <- as.matrix(py_to_r(counts))

sample_names <- fread(args[3], sep="\n", header=FALSE)

# Ensure that sample_names is a vector, not a data table
sample_names <- sample_names$V1

storage.mode(counts) <- "integer"
colnames(counts) <- sample_names

# metadata <- read_delim(args[4], delim='\t', col_names=T)
# #rownames(metadata) <- metadata$id
colData <- data.frame(row.names=sample_names, sample_ids=sample_names)

print('Making DESeq dataset')
dds <- DESeqDataSet(SummarizedExperiment(assays = list(counts=counts), colData=colData), design=~1)
# dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~1)
rm(counts)
gc()

# Provide NULL or non-existent norm_factors file for conventional VST
if (is.null(args[2]) | file.exists(args[2])) {
  print('Reading norm factors')
  norm_factors <- as.matrix(py_to_r(np$load(args[2])))
  print('Converting norm factors to R')
  normalizationFactors(dds) <- norm_factors
  rm(norm_factors)
  gc()
} else {
  print("Calculating size factors")
  dds <- estimateSizeFactors(dds)
}

if (is.null(params_f)) {
    print('Calculating and saving VST params')
    nsub <- 1000
    fitType <- "parametric"
    # code below was copied from https://github.com/mikelove/DESeq2/blob/master/R/vst.R
    # dispersionFunction is not getting saved otherwise
    baseMean <- MatrixGenerics::rowMeans(counts(object, normalized=TRUE))
    if (sum(baseMean > 5) < nsub) {
        stop("less than 1000 rows with mean normalized count > 5, 
        it is recommended to use varianceStabilizingTransformation directly")
    }

    # subset to a specified number of genes with mean normalized count > 5
    object.sub <- object[baseMean > 5,]
    baseMean <- baseMean[baseMean > 5]
    o <- order(baseMean)
    idx <- o[round(seq(from=1, to=length(o), length=nsub))]
    object.sub <- object.sub[idx,]

    # estimate dispersion trend
    object.sub <- estimateDispersionsGeneEst(object.sub, quiet=TRUE)
    object.sub <- estimateDispersionsFit(object.sub, fitType=fitType, quiet=TRUE)

    # assign to the full object
    dispersionFunction(object) <- dispersionFunction(object.sub)
} else {
    print('Use existing VST params')
    df <- readRDS(params_f)
    dispersionFunction(dds) <- df
}

params_file_name <- paste(prefix, ".params.RDS", sep='')
if (file.exists(params_file_name)) {
  print(paste('Parameters were not saved. File ', params_file_name, ' exists.', sep=''))
} else {
  saveRDS(dispersionFunction(dds), file=params_file_name)
}
# vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
# np$save(args[2], np$array(assay(vsd), dtype='float32'))

dds_name <- paste(prefix, ".dds.RDS", sep='')
timing <- system.time({
  saveRDS(dds, dds_name, compress = FALSE)
})

print(timing)