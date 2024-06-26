library(variancePartition)
library(data.table)
library(rhdf5)


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


data <- h5read(file_path, 'vst', 
    start=c(1, start_index), 
    count=c(nrow(meta), count)
)
data <- t(data)

sample_names <- h5read(file_path, 'sample_names')
colnames(data) <- sample_names
row.names(data) <- row.names(dhs_meta)

# Sort the DataFrame according to the sample_names
sorted_metadata <- meta[sample_names, ]

formula <- args[6]


print('Fitting model')

varPart <- fitExtractVarPartModel(data, formula, sorted_metadata)

write.table(varPart , "tmp.txt", sep="\t", row.names=FALSE, quote = FALSE)
vp <- fread("tmp.txt")
stopifnot(all(identical(row.names(varPart), row.names(dhs_meta))))
write.table(cbind(dhs_meta, vp) , args[7], sep="\t", row.names=FALSE, quote = FALSE)
