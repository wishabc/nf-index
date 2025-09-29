library(variancePartition)
library(data.table)
library(reticulate)
np <- import("numpy", convert=FALSE)

print("Finished imports")
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
  stop("At least 8 input arguments should be supplied", call.=FALSE)
}

sample_names <- fread(args[1], sep="\n", header=FALSE)$V1

sample_meta <- as.data.frame(fread(args[2], stringsAsFactors = TRUE))
row.names(sample_meta) <- sample_meta$ag_id
sample_meta <- sample_meta[sample_names, ]


start_index <- as.integer(args[3])
count <- as.integer(args[4])
matrix_path <- args[5]
dhs_meta <- fread(args[6])

count <- min(count, nrow(dhs_meta) - start_index + 1)

dhs_meta <- dhs_meta[start_index:(start_index + count - 1), ]
row.names(dhs_meta) <- dhs_meta$chunk_id

print(start_index)
print(start_index - 1 + count)
np_array <- np$load(matrix_path, mmap_mode = 'r')[(start_index - 1):(start_index - 1 + count)]

print("Read as np")
data <- py_to_r(np_array)

print("Data loaded")
print(dim(data))
print(dim(sample_meta))
print(dim(dhs_meta))
colnames(data) <- row.names(sample_meta)
row.names(data) <- row.names(dhs_meta)

formula <- args[7]

print("Fitting models")
good <- apply(data, 1, function(x) var(x)) > 0


varPart_sub <- fitExtractVarPartModel(data[good, ], formula, sample_meta)
varPart <- varPart_sub[FALSE, ]
varPart <- varPart[rep(NA, nrow(data)), ]
rownames(varPart) <- rownames(data)
varPart[good, ] <- varPart_sub

vp_df <- as.data.frame(varPart)
stopifnot(all(identical(row.names(vp_df), row.names(dhs_meta))))
write.table(cbind(dhs_meta, vp_df), args[8], sep="\t", row.names=FALSE, quote = FALSE)
