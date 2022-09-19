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

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("At least four input files should be supplied", call.=FALSE)
}

counts <- np$load(args[1])
norm_factors <- np$load(args[2])
sample_names <- fread(args[3], header=FALSE)
sample_names <- data.table(t(sapply(sample_names, FUN=function(x) { sprintf('AG_%s', x)})))

counts <- as.data.frame(counts, stringsAsFactors = F)
colnames(counts) <- sample_names

metadata <- read_delim(metadata, delim = '\t', col_names=T)
rownames(metadata) <- metadata$indiv_id

# Applying DESEQ with norm_factors
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
normalizationFactors(dds) <- norm_factors
vsd <- vst(dds, blind = F)
np$save(paste(args[5], ".vst.txt", sep=''), np$array(assay(vsd)))

# Applying DESEQ no norm_factors
dds2 <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
vsd2 <- vst(dds2, blind = F)

np$save(paste(args[5], ".no_sf.vst.txt", sep=''), np$array(assay(vsd2)))
