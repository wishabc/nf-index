library(variancePartition)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("At least three input arguments should be supplied", call.=FALSE)
}

meta = read.delim(args[1])
rownames(meta) <- meta$ag_id
meta$sample_id <- meta$ag_id
meta <- as.data.frame(meta)

data <- read.delim(args[2], header=FALSE)
data <- as.data.frame(data, stringsAsFactors = F)
sample_names <- fread(args[3], sep="\n", header=FALSE)
sample_names <- sample_names$V1

colnames(data) <- sample_names


form <- ~ dedupped_subsampled_spot1 + log(read_depth) # + (1|sex) + 

varPart <- fitExtractVarPartModel(data, form, meta)
vp <- sortCols(varPart)
write.delim(vp, args[4])


