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

vsd <- varianceStabilizingTransformation(dds, blind = F)
np$save(args[2], np$array(vsd, dtype='float32'))
