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
sample_names <- data.table(sample_names)

counts <- as.data.frame(counts, stringsAsFactors = F)
colnames(counts) <- sample_names

metadata <- read_delim(args[4], delim = '\t', col_names=T)
rownames(metadata) <- metadata$uniq_id
print(counts)
print(metadata)
# Applying DESEQ with norm_factors
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
normalizationFactors(dds) <- norm_factors
vsd <- vst(dds, blind = F)
rm(dds)
gc()
np$save(paste(args[5], ".vst.txt", sep=''), np$array(assay(vsd)))
rm(vsd)
gc()
# Applying DESEQ no norm_factors
dds2 <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~1)
rm(counts)
gc()
vsd2 <- vst(dds2, blind = F)
rm(dds2)
gc()
np$save(paste(args[5], ".no_sf.vst.txt", sep=''), np$array(assay(vsd2)))
