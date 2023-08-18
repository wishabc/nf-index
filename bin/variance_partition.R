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

meta <- meta[match(sample_names, row.names(meta)), ]


formula <- ~ dedupped_subsampled_spot1 + log(read_depth) + dupRate_5M + (1 | donor_sex) + (1 | library_kit) + (1 | short_ontology)

print('Fitting model')
safeFitExtractVarPartModel <- function(data, formula, meta) {
  # Create an empty list to store the results for each row
  result_list <- vector("list", nrow(data))
  num_fixed_effects <- length(all.vars(formula))
  num_random_effects <- length(attr(terms(formula, random = ~ 1), "specials")$re)
  total_NA <- num_fixed_effects + num_random_effects + 1
  
  # Iterate through the rows (genes)
  for (i in 1:nrow(data)) {
    # Extract the row
    cat("Processing row:", i, "\n")
    row_data <- data[i, , drop=FALSE]
    # Try fitting the model to the row, and handle any errors
    result_list[[i]] <- tryCatch(
      {
        varPart <- fitExtractVarPartModel(row_data, formula, meta)
        if (isSingular(varPart)) {
            warning(paste("Singular fit encountered in row", i))
            varPart <- NA # or another appropriate value
        }
        cat("Row processed successfully:", i, "\n")
        return(varPart)
      },
      error = function(e) {

        warning(paste("Singular fit error encountered in row", i, ":", e))
        na_result <- data.frame(t(rep(NA, total_NA))) # Return NA for each element in the formula + one for residuals
        rownames(na_result) <- rownames(row_data) # Set the row names of the NA result to match row_data
        cat("Returning NA result for row:", i, "\n")
        return(na_result)
      }
    )
  }
  cat("Number of error rows:", length(error_rows), "\n") # Print the number of error rows

  # Combine the results into a final object (you may need to adjust this part based on the desired format)
  final_result <- do.call(rbind, result_list)
  return(final_result)
}

varPart <- safeFitExtractVarPartModel(data, formula, meta)

stopifnot(identical(row.names(varPart), row.names(dhs_meta)))
write.table(cbind(dhs_meta, varPart) , args[6], sep="\t", row.names=FALSE, quote = FALSE)


