import numpy as np 
import os
import sys
import pdb







# Get mapping from clusterID to array of genes mapped to that cluster
def get_cluster_to_gene_mapping(cluster_info_file):
	# Initialize mapping from cluster to gene array
	mapping = {}
	# Used to skip header
	head_count = 0
	# Stream file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relvent fields from line
		cluster_id = data[0]
		gene_string = data[2]
		# Get array of genes from string
		gene_arr = gene_string.split(',')
		# Add cluster to mapping dictionary (if it has never been seen before)
		if cluster_id not in mapping:
			mapping[cluster_id] = []
		else:
			print('twice!!')
			pdb.set_trace()
		# Add genes
		for gene in gene_arr:
			mapping[cluster_id].append(gene)
	f.close()
	# Mapping sure each cluster only maps to a UNIQUE set of genes
	for cluster_id in mapping.keys():
		mapping[cluster_id] = np.unique(mapping[cluster_id])
	return mapping

# Make data structure where each key is a gene and maps to:
# a. An integer that counts number of clusters associated with that gene
# b. A vewctor of length number of samples (in this tissue) where each element is the min(pvalue) across all clusters for the given gene
def get_gene_based_pvalue_data_structure(cluster_level_outlier_file, cluster_to_gene_mapping):
	# Initialize data structure
	gene_level_data_structure = {}
	# Used to skip header
	head_count = 0
	# Stream cluster_level_outlier_file
	f = open(cluster_level_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# skip header
		if head_count == 0:
			head_count = head_count + 1
			# Count number of samples
			sample_names = np.asarray(data[1:])
			num_samples = len(sample_names)
			continue
		# Standard line
		# extract relevent fields
		cluster_id = data[0]
		pvalue_vector = np.asarray(data[1:]).astype(float)
		# Get array of genes that this cluster maps to
		gene_arr = cluster_to_gene_mapping[cluster_id]
		# Loop through all genes
		for gene in gene_arr:
			# Add gene to data structure if it has never been seen before
			if gene not in gene_level_data_structure:
				gene_level_data_structure[gene] = {}
				gene_level_data_structure[gene]['num_clusters'] = 0  # a. An integer that counts number of clusters associated with that gene
				gene_level_data_structure[gene]['pvalue'] = np.ones(num_samples)  # b. A vewctor of length number of samples (in this tissue) where each element is the min(pvalue) across all clusters for the given gene
			# Update gene level data structure for current cluster
			gene_level_data_structure[gene]['num_clusters'] = gene_level_data_structure[gene]['num_clusters'] + 1
			gene_level_data_structure[gene]['pvalue'] = np.minimum(gene_level_data_structure[gene]['pvalue'], pvalue_vector)
	f.close()
	return gene_level_data_structure, sample_names

# Correct pvalues (accounting for the number of clusters we are taking the minimum over)
def correct_pvalues(uncorrected_pvalues, num_clusters):
	# Initialize output vector
	corrected_pvalues = []
	# Loop through samples
	for uncorrected_pvalue in uncorrected_pvalues:
		# Correct pvalues for min transformation
		corrected_pvalue = 1.0 - ((1.0 - uncorrected_pvalue)**num_clusters)
		# Add to array
		corrected_pvalues.append(corrected_pvalue)
	return np.asarray(corrected_pvalues)

# For each gene, correct the gene level pvalues (accounting for the number of clusters we are taking the minimum over)
# Then print to output file
def correct_gene_level_pvalues_and_print_to_output_file(gene_level_outlier_file, gene_based_data_structure, sample_names):
	# Open output file handle
	t = open(gene_level_outlier_file, 'w')
	# Print header
	t.write('GENE_ID\t' + '\t'.join(sample_names) + '\n')
	# Loop through genes
	for gene in gene_based_data_structure.keys():
		# Extract information from gene based data structure
		uncorrected_pvalues = gene_based_data_structure[gene]['pvalue']
		num_clusters = gene_based_data_structure[gene]['num_clusters']
		# Correct pvalues (accounting for the number of clusters we are taking the minimum over)
		corrected_pvalues = correct_pvalues(uncorrected_pvalues, num_clusters)
		# Print to output file
		t.write(gene + '\t' + '\t'.join(corrected_pvalues.astype(str)) + '\n')
	# Close output file
	t.close()


###########################
# Command line arguments
###########################
output_root = sys.argv[1]  # Input and output root
cluster_info_file = sys.argv[2]  # File containing junction/cluster/gene mapping information



#  Input file containing oultier calls at the cluster level (for a specific tissue)
cluster_level_outlier_file = output_root + 'merged_emperical_pvalue.txt'

#  Output file containing oultier calls at the gene level (for a specific tissue)
gene_level_outlier_file = output_root + 'merged_emperical_pvalue_gene_level.txt'





# Get mapping from clusterID to array of genes mapped to that cluster
cluster_to_gene_mapping = get_cluster_to_gene_mapping(cluster_info_file)


# Make data structure where each key is a gene and maps to:
# a. An integer that counts number of clusters associated with that gene
# b. A vewctor of length number of samples (in this tissue) where each element is the min(pvalue) across all clusters for the given gene
# Also extract vector of sample names
gene_based_data_structure, sample_names = get_gene_based_pvalue_data_structure(cluster_level_outlier_file, cluster_to_gene_mapping)


# For each gene, correct the gene level pvalues (accounting for the number of clusters we are taking the minimum over)
# Then print to output file
correct_gene_level_pvalues_and_print_to_output_file(gene_level_outlier_file, gene_based_data_structure, sample_names)


