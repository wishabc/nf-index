library(readr)
library(dplyr)
library(reticulate)
library(stringr)
library(data.table)
library(DESeq2)
library(SummarizedExperiment)
library(tibble)
library(purrr)
library(ggplot2)

contrast_species_vs_mean <- function(dds, ann, vs_mean = TRUE) {
  coefs   <- resultsNames(dds)
  coefs_h <- grep("\\.human$", coefs, value = TRUE)
  coefs_m <- grep("\\.mouse$", coefs, value = TRUE)
  anns    <- sub("^ann_species", "", sub("\\.(human|mouse)$", "", coefs_h))
  J       <- length(anns)

  if (!ann %in% anns) {
    stop("`ann` must be one of: ", paste(anns, collapse = ", "))
  }

  # start with zeros
  w <- setNames(numeric(length(coefs)), coefs)

  # optionally subtract/add the mean species effect
  if (vs_mean && J > 0) {
    w[coefs_h] <- -1 / J
    w[coefs_m] <-  1 / J
  }
  # add back the chosen annotation's (human - mouse)
  h <- paste0("ann_species", ann, ".human")
  m <- paste0("ann_species", ann, ".mouse")
  w[h] <- w[h] + 1
  w[m] <- w[m] - 1

  w
}

export_deseq_data <- function(dds, prefix) {
  stopifnot(inherits(dds, "DESeqDataSet"))

  ## 2) global interaction LRT (reuse if already run)
  print('global LRT')
  dds <- nbinomLRT(
    dds,
    full    = ~ 0 + extended_annotation * species,
    reduced = ~ 0 + extended_annotation + species
  )
  res_lrt <- results(dds)
  lrt_tbl <- data.frame(
    gene     = rownames(dds),
    baseMean = res_lrt$baseMean,
    stat     = res_lrt$stat,
    pvalue   = res_lrt$pvalue,
    padj     = res_lrt$padj
  )
  write.table(lrt_tbl, paste0("LRT.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (!("mu" %in% assayNames(dds))) stop("Assay 'mu' missing. Run nbinomWaldTest() or nbinomLRT() first.")
    

  ## 1) per-sample long data (box/strip + fitted means)
  print('Plot data')
  nc   <- (counts(dds, normalized = FALSE) + 1) / normalizationFactors(dds)
  mu   <- assay(dds, "mu")
  q <- (mu + 1) / normalizationFactors(dds)
  gids <- rownames(dds); sids <- colnames(dds)
  stopifnot(ncol(dds) > 0)

  df_plot <- data.frame(
    gene                 = rep(gids, times = length(sids)),
    sample               = rep(sids, each  = length(gids)),
    log10_norm_count_p1  = log10(as.numeric(nc)),
    log10_q_hat_p1      = log10(as.numeric(q))
  )
  md <- as.data.frame(colData(dds))
  keep_cols <- intersect(c("species","extended_annotation"), colnames(md))
  md$sample <- rownames(md)
  df_plot <- merge(df_plot, md[, c(keep_cols, "sample")], by = "sample", all.x = TRUE, sort = FALSE)
  write.table(df_plot, paste0("plot_data.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    

  ## 3) ALL Wald coefficients specie*extended_annotation vs 0
  print('Wald raw betas')
  design(dds) <- ~ 0 + ann_species
  dds_w <- nbinomWaldTest(dds)
  rn <- resultsNames(dds_w)

  all_coef <- do.call(rbind, lapply(rn, function(nm) {
    r <- results(dds_w, name = nm)
    data.frame(
      gene           = rownames(dds_w),
      coefficient    = nm,
      log2FoldChange = r$log2FoldChange,
      lfcSE          = r$lfcSE,
      stat           = r$stat,
      pvalue         = r$pvalue,
      padj           = r$padj,
      stringsAsFactors = FALSE
    )
  }))
  fn_all <- paste0("Wald_coeffs_ALL.", prefix, ".tsv")
  write.table(all_coef, fn_all, sep = "\t", quote = FALSE, row.names = FALSE)

  ## 4) Interaction significance species1 - species2 for a given cell type vs mean across cell types
  print('Wald specie diff for each CT vs mean')
  anns <- sub("^ann_species(.*)\\.human$", "\\1", grep("\\.human$", resultsNames(dds_w), value=TRUE))
  all_coef_contrast <- do.call(rbind, lapply(anns, function(nm) {
    r <- results(dds_w, contrast = contrast_species_vs_mean(dds_w, nm))
    data.frame(
      gene           = rownames(dds_w),
      annotation     = nm,
      log2FoldChange = r$log2FoldChange,
      lfcSE          = r$lfcSE,
      stat           = r$stat,
      pvalue         = r$pvalue,
      padj           = r$padj,
      stringsAsFactors = FALSE
    )
  }))
  write.table(all_coef_contrast, paste0("Wald_coeffs_INTERACTION.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ## 4.5) Interaction species1 - species2 for a given cell type vs 0
  print('Wald specie diff for each CT vs 0')
  anns <- sub("^ann_species(.*)\\.human$", "\\1", grep("\\.human$", resultsNames(dds_w), value=TRUE))
  all_coef_contrast_zero <- do.call(rbind, lapply(anns, function(nm) {
    r <- results(dds_w, contrast = contrast_species_vs_mean(dds_w, nm, vs_mean=FALSE))
    data.frame(
      gene           = rownames(dds_w),
      annotation     = nm,
      log2FoldChange = r$log2FoldChange,
      lfcSE          = r$lfcSE,
      stat           = r$stat,
      pvalue         = r$pvalue,
      padj           = r$padj,
      stringsAsFactors = FALSE
    )
  }))
  write.table(all_coef_contrast_zero, paste0("Wald_coeffs_INTERACTION_ZERO.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    
  ## 5) species1 - species2 average
  print('Wald specie avg')
  coefs <- resultsNames(dds_w)
  coefs_h <- grep("\\.human$", coefs, value=TRUE)
  coefs_m <- grep("\\.mouse$", coefs, value=TRUE)
  w <- setNames(rep(0, length(coefs)), coefs)
  w[coefs_h] <-  1 / length(coefs_h)
  w[coefs_m] <- -1 / length(coefs_m)
  r <- results(dds_w, contrast = w)
  r <- data.frame(
      gene           = rownames(dds_w),
      log2FoldChange = r$log2FoldChange,
      lfcSE          = r$lfcSE,
      stat           = r$stat,
      pvalue         = r$pvalue,
      padj           = r$padj,
      stringsAsFactors = FALSE
    )
  write.table(r, paste0("Wald_coeffs_SPECIES.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
  ## 6) Sample data
  print('Col data')
  write.table(colData(dds), paste0("samples_data.", prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(plot = df_plot, LRT = lrt_tbl, Wald_ALL = all_coef))
}

get_chunk <- function(dds, total_chunks, i) {
  n_genes <- nrow(dds)
  chunk_size <- ceiling(n_genes / total_chunks)

  start <- (i - 1) * chunk_size + 1
  end   <- min(i * chunk_size, n_genes)

  dds[start:end, ]
}

args = commandArgs(trailingOnly=TRUE)

print('Loading dds')
dds_chunk <- readRDS(args[2])
print('Slicing chunk')
dds_chunk <- get_chunk(dds_chunk, total_chunks = as.integer(args[3]), i = as.integer(args[4]))
print('Estimating point dispersions')
dds_chunk <- estimateDispersionsGeneEst(dds_chunk)
print('Estimating dispersions MAP')
dds_chunk <- estimateDispersionsMAP(dds_chunk)
print('Exporting DESeq data')
export_deseq_data(dds_chunk, prefix = args[1])
