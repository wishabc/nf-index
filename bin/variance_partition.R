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

data <- t(
        h5read(args[4], 
        'vst', 
        start=c(1, start_index), 
        count=c(ncol(meta), count)
    )
)

sample_names <- h5read(args[4], 'sample_names')


colnames(data) <- sample_names


form <- ~ dedupped_subsampled_spot1 + log(read_depth) # + (1|sex) + 

varPart <- fitExtractVarPartModel(data, form, meta)
vp <- sortCols(varPart)
write.delim(vp, args[4])


