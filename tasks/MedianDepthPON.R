#!/usr/bin/Rscript

library(data.table)	

# Read input file paths from command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Initialize a list to store the matrices
matrices <- list()

# Read each file and store it in the list
for (i in 1:length(args)) {
  matrices[[i]] <- as.matrix(read.table(args[i]))
}

# Combine matrices into an array
all_matrices <- array(unlist(matrices), dim = c(nrow(matrices[[1]]), ncol(matrices[[1]]), length(matrices)))

# Compute median matrix
median_matrix <- apply(all_matrices, c(1, 2), median)

# save the median matrix
write.table(all_matrices, file = "all_matrices_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(median_matrix, file = "PON_median_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
