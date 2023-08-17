library(variancePartition)
library(data.table)
library(rhdf5)


args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least three input arguments should be supplied", call.=FALSE)
}

meta = read.delim(args[1])
rownames(meta) <- meta$ag_id
meta$sample_id <- meta$ag_id
meta <- as.data.frame(meta)

start_index <- as.integer(args[2])
count <- as.integer(args[3])
file_path <- args[4]

info <- h5ls(file_path)
vst_info <- info[info$name == "vst", ]
num_dhs <- as.numeric(unlist(strsplit(vst_info$dim, " x ")))[2]


data <- h5read(file_path, 'vst', 
    start=c(1, start_index), 
    count=c(nrow(meta), min(count, num_dhs - start_index + 1))
)
data <- t(data)

sample_names <- h5read(file_path, 'sample_names')
colnames(data) <- sample_names


meta <- meta[match(sample_names, rownames(meta)), ]




form <- ~ dedupped_subsampled_spot1 + log(read_depth) # + (1|sex) + 

varPart <- fitExtractVarPartModel(data, form, meta)
write.table(varPart , args[5], sep="\t", row.names=FALSE, quote = FALSE)


