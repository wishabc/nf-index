library(variancePartition)
library(data.table)
library(reticulate)
np <- import("numpy", convert=FALSE)

print("Finished imports")
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 7) {
  stop("At least six input arguments should be supplied", call.=FALSE)
}

meta <- read.delim(args[1])
row.names(meta) <- meta$ag_id
meta$sample_id <- meta$ag_id
meta <- as.data.frame(meta)

start_index <- as.integer(args[2])
count <- as.integer(args[3])
file_path <- args[4]
dhs_meta <- fread(args[5])

count <- min(count, nrow(dhs_meta) - start_index + 1)

dhs_meta <- dhs_meta[start_index:(start_index + count - 1), ]
row.names(dhs_meta) <- dhs_meta$chunk_id


data <- np$load(file_path, mmap_mode = 'r')[start_index: (start_index + count - 1), ]

print("Data loaded")
print(dim(data))
print(dim(meta))
print(dim(dhs_meta))
colnames(data) <- row.names(meta)
row.names(data) <- row.names(dhs_meta)

formula <- args[6]


print("Fitting models")

varPart <- fitExtractVarPartModel(data, formula, meta)

vp_df <- as.data.frame(varPart)
stopifnot(all(identical(row.names(vp_df), row.names(dhs_meta))))
write.table(cbind(dhs_meta, vp_df), args[7], sep="\t", row.names=FALSE, quote = FALSE)
