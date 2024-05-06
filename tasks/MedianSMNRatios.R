#!/usr/bin/Rscript

library(data.table)	

# Read input file paths from command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Initialize a list to store the matrices
matrices <- list()

# Read each file and store it in the list
for (i in 1:length(args)) {
  matrices[[i]] <- as.matrix(read.table(args[i]), header=FALSE)
}

# Combine matrices into an array. ncol = #PONs x #norm sites
all_matrices <- array(unlist(matrices), dim = c(nrow(matrices[[1]]), ncol(matrices[[1]]), length(matrices)))

# Compute median matrix of ratios per SNP. resulting median matrix: ncol = #norm sites
median_matrix_perSNP <- apply(all_matrices, c(1, 2), median)

# Compute median matrix of ratios across all SNPs
median_ratios <- apply(median_matrix_perSNP, 1, median)

# save the median matrix
write.table(median_matrix_perSNP, file = "median_matrix_perSNP.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(median_ratios, file = "median_ratios.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
