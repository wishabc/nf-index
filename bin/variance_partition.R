library(variancePartition)
library(data.table)
library(rhdf5)


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("At least six input arguments should be supplied", call.=FALSE)
}

meta = read.delim(args[1])
rownames(meta) <- meta$ag_id
meta$sample_id <- meta$ag_id
meta <- as.data.frame(meta)

start_index <- as.integer(args[2])
count <- as.integer(args[3])
file_path <- args[4]
dhs_meta <- read.delim(args[5])

end_index <- min(count, nrow(dhs_meta) - start_index + 1)

dhs_meta <- dhs_meta[start_index:end_index, ]
row.names(dhs_meta) <- dhs_meta$chunk_id


data <- h5read(file_path, 'vst', 
    start=c(1, start_index), 
    count=c(nrow(meta), end_index)
)
data <- t(data)

sample_names <- h5read(file_path, 'sample_names')
colnames(data) <- sample_names
print(dim(data), dim(dhs_meta))
row.names(data) <- row.names(dhs_meta)

meta <- meta[match(sample_names, row.names(meta)), ]


form <- ~ dedupped_subsampled_spot1 + log(read_depth) + (1|donor_sex)
print('Fitting model')
varPart <- fitExtractVarPartModel(data, form, meta)
stopifnot(identical(row.names(varPart), row.names(dhs_meta)))
write.table(varPart , args[6], sep="\t", row.names=FALSE, quote = FALSE)


