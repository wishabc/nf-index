library(readr)
library(dplyr)

library(stringr)
library(data.table)
library(DESeq2)
library(SummarizedExperiment)
library(BiocParallel)

library(reticulate)
ad <- import("anndata")

args = commandArgs(trailingOnly=TRUE)


anndata <- ad$read_zarr(args[1])
anndata

sample_meta <- anndata$obs
dhs_meta <- anndata$var
counts <- t(anndata$layers['counts'])
storage.mode(counts) <- "integer"
row.names(counts) <- row.names(dhs_meta)
colnames(counts) <- row.names(sample_meta)


dds <- DESeqDataSet(
    SummarizedExperiment(assays=list(counts=counts), colData=sample_meta),
    design=~species
)
print('DEseq dataset done')
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, parallel=FALSE, minReplicatesForReplace=5)

res_df <- as.data.frame(results(dds))
res_df$dhs_id <- rownames(res_df)
res_df$prefix <- args[2]

write.table(
    res_df,
    file=args[3],
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)
