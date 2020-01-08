args = commandArgs(trailingOnly=TRUE)




compute_tpm <- function(raw_counts, gene_length) {
	x <- raw_counts / gene_length
	tpm_mat <- t( t(x) * 1e6 / colSums(x) )
	return(tpm_mat)
}

# Helper method to save DGE data structure to tab-deliminated text file
save_matrix <- function(counts, output_file) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0("Sample_name\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names,sep="\t")

}

raw_counts_file <- args[1]  # input file
count_file <- args[2] # output file
tpm_file <- args[3]  # output file
log_tpm_file <- args[4]  # output file



# Load in raw counts data
raw_counts_data <- read.table(raw_counts_file, header=TRUE)

# Matrix of just raw counts
raw_counts = raw_counts_data[,7:(dim(raw_counts_data)[2])]

# length of genes
gene_length = raw_counts_data[,6]

# Names of genes
gene_ids = raw_counts_data[,1]

colnames(raw_counts) = sub("X", "", colnames(raw_counts))
# Compute TPM mat
tpm_mat <- compute_tpm(raw_counts, gene_length)

# Save tpm mat
rownames(raw_counts) = gene_ids
save_matrix(raw_counts, count_file)

# Save tpm mat
rownames(tpm_mat) = gene_ids
save_matrix(tpm_mat, tpm_file)

# Save log tpm mat
log_tpm_mat <- log2(tpm_mat + 2)
save_matrix(log_tpm_mat, log_tpm_file)
