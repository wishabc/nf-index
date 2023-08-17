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

hdf5_file <- H5Fopen(args[4])
data <- h5read(hdf5_file, 
    'vst', 
    start=c(start_index, 1), 
    count=c(count, -1)
)

sample_names <- h5read(hdf5_file, 'sample_names')

H5Fclose(hdf5_file)

colnames(data) <- sample_names


form <- ~ dedupped_subsampled_spot1 + log(read_depth) # + (1|sex) + 

varPart <- fitExtractVarPartModel(data, form, meta)
vp <- sortCols(varPart)
write.delim(vp, args[4])


