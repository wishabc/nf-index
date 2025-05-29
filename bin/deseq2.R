#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(reticulate)
library(stringr)
library(data.table)
library(DESeq2)

np <- import("numpy")

args = commandArgs(trailingOnly=TRUE)

dds <- readRDS(args[1])
n_genes <- nrow(dds)
n_samples <- ncol(dds)
chunk_size <- 100000
# Initialize empty matrix for transformed results
vsd <- matrix(NA_real_, nrow=n_genes, ncol=n_samples)

# Process in row chunks
for (start in seq(1, n_genes, by=chunk_size)) {
  end <- min(start + chunk_size - 1, n_genes)
  cat("Processing DHSs", start, "to", end, "\n")
  
  dds_chunk <- dds[start:end, ]
  vsd_chunk <- varianceStabilizingTransformation(dds_chunk, blind = FALSE)
  vsd[start:end, ] <- assay(vsd_chunk)
}
#vsd <- varianceStabilizingTransformation(dds, blind = F)
np$save(args[2], np$array(vsd, dtype='float32'))
