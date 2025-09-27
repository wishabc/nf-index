#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(reticulate)
library(stringr)
library(data.table)
library(DESeq2)

np <- import("numpy")

args = commandArgs(trailingOnly=TRUE)


prefix <- args[1]

formula <- as.formula(args[3])

dds <- readRDS(args[2])
design(dds) <- formula


if (length(args) >= 4) {
    print("Taking existing params")
    params_f <- args[4]
    print('Use existing VST params')
    df <- readRDS(params_f)
    dispersionFunction(dds) <- df
} else {
    print('Calculating and saving VST params')
    nsub <- 1000
    min_baseMean <- 2
    fitType <- "parametric"
    object <- dds

    # code below was copied from https://github.com/mikelove/DESeq2/blob/master/R/vst.R
    # dispersionFunction is not getting saved otherwise
    baseMean <- rowMeans(counts(object, normalized=TRUE))
    if (sum(baseMean > min_baseMean) < nsub) {
        stop("less than 1000 rows with mean normalized count > 5, 
        it is recommended to use varianceStabilizingTransformation directly")
    }

    # subset to a specified number of genes with mean normalized count > 5
    object.sub <- object[baseMean > min_baseMean,]
    baseMean <- baseMean[baseMean > min_baseMean]
    o <- order(baseMean)
    idx <- o[round(seq(from=1, to=length(o), length=nsub))]
    object.sub <- object.sub[idx,]

    # estimate dispersion trend
    object.sub <- estimateDispersionsGeneEst(object.sub, quiet=TRUE)
    object.sub <- estimateDispersionsFit(object.sub, fitType=fitType, quiet=TRUE)

    object.sub <- estimateDispersionsMAP(object.sub)

    dispersion_plot_name <- paste(prefix, ".dispersions_plot.pdf", sep='')
    pdf(dispersion_plot_name, width = 3, height = 3)
    plotDispEsts(object.sub, ylim = c(1e-4, 10), xlim = c(2, 2000))
    dev.off()

    # assign to the full object
    dispersionFunction(object) <- dispersionFunction(object.sub)
}

params_file_name <- paste(prefix, ".dispersion_function.RDS", sep='')
if (file.exists(params_file_name)) {
  print(paste('Parameters were not saved. File ', params_file_name, ' exists.', sep=''))
} else {
  saveRDS(dispersionFunction(dds), file=params_file_name)
}

vsd <- varianceStabilizingTransformation(dds, blind = F)
vsd_name <- paste(prefix, ".vst.npy", sep='')
np$save(vsd_name, np$array(vsd, dtype='float32'))
