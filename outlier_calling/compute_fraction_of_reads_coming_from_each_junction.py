import numpy as np 
import os
import sys
import pdb




def extract_raw_cluster_jxn_data_structure(tissue_specific_jxn_file):
	# Used to skip header
	head_count = 0
	# Initialize cluster_jxn_data_structure
	cluster_jxn_data_structure = {}
	# Stream input file
	f = open(tissue_specific_jxn_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			samples = np.asarray(data[1:])
			continue
		# Standard line
		# Extract Jxn info from current line
		key = data[0]
		key_info = key.split(':')
		chrom_num = key_info[0]
		start = key_info[1]  # 5' ss
		end = key_info[2]  # 3' ss
		cluster_id = key_info[3]  # Name of cluster
		genes = key_info[4].split(',')  # Array of genes that cluster is mapped to
		jxn_read_counts = np.asarray(data[1:]).astype(float)  # Vector of raw read counts for this junction

		# Add jxn to  cluster_jxn_data_structure
		if cluster_id not in cluster_jxn_data_structure:  # If we've never seen this cluster before
			cluster_jxn_data_structure[cluster_id] = []
		cluster_jxn_data_structure[cluster_id].append(jxn_read_counts)  # Add read counts from current junction
	f.close()
	# Convert from list of arrays to matrix (in each cluster)
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		# Create matrix from list of arrays
		jxn_matrix = np.transpose(np.asmatrix(cluster_jxn_data_structure[cluster_id]))
		# Add this new matrix to the data structure
		cluster_jxn_data_structure[cluster_id] = {}
		cluster_jxn_data_structure[cluster_id]['samples'] = samples
		cluster_jxn_data_structure[cluster_id]['jxn_matrix'] = jxn_matrix
	return cluster_jxn_data_structure, samples

# Remove clusters with more than $max_number_of_junctions_per_cluster
def max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_number_of_junctions_per_cluster):
	# Initialize new cluster_jxn_data_structure
	new_jxn_structure = {}
	# Loop through clusters
	for cluster_id in cluster_jxn_data_structure.keys():
		samples = cluster_jxn_data_structure[cluster_id]['samples']
		jxn_mat = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# Get dimensionality of this cluster
		N, K = jxn_mat.shape
		# If there are more than max_number_of_junctions_per_cluster
		if K > max_number_of_junctions_per_cluster:
			pass
		else: # Less than or equal
			new_jxn_structure[cluster_id] = {}
			new_jxn_structure[cluster_id]['samples'] = samples
			new_jxn_structure[cluster_id]['jxn_matrix'] = jxn_mat
	return new_jxn_structure

# Remove samples with less than $min_reads_per_sample_in_cluster reads (summed acrosss all junctions in cluster)
def min_reads_per_sample_filter(cluster_jxn_data_structure, min_reads_per_sample_in_cluster):
	# Initialize new data structure (post filtering)
	new_jxn_structure = {}
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		# Extract relevent fields from cluster
		samples = cluster_jxn_data_structure[cluster_id]['samples']
		jxn_mat = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# get dimensionality of cluster
		N, K = jxn_mat.shape
		# Itialize vectors of samples and jxn mat (post filtering)
		new_samples = []
		new_jxn_mat = []
		# Loop through samples
		for n in range(N):
			# Vector spanning this samples junction counts for this cluster
			sample_counts = np.squeeze(np.asarray(jxn_mat[n,:]))
			sample_id = samples[n]
			# Sum of reads across all junctions for this cluster
			total_sample_read_counts = np.sum(sample_counts)
			# Check if sample passes filter
			if total_sample_read_counts >= min_reads_per_sample_in_cluster:
				new_samples.append(sample_id)
				new_jxn_mat.append(sample_counts)
		# Convert from list of arrays to matrices
		new_jxn_mat = np.asmatrix(new_jxn_mat)
		new_samples = np.asarray(new_samples)
		# Add filtered matrices to jxn_structure
		new_jxn_structure[cluster_id] = {}
		new_jxn_structure[cluster_id]['samples'] = new_samples
		new_jxn_structure[cluster_id]['jxn_matrix'] = new_jxn_mat
	return new_jxn_structure


def add_intercept_covariate_to_cluster_jxn_data_structure(cluster_jxn_data_structure):
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		num_samples = len(cluster_jxn_data_structure[cluster_id]['samples'])
		cov_mat = np.ones((num_samples, 1))
		cluster_jxn_data_structure[cluster_id]['covariate_matrix'] = cov_mat
	return cluster_jxn_data_structure


# Add Covariate Matrix to cluster_jxn_data_structure
def add_covariates_to_cluster_jxn_data_structure(cluster_jxn_data_structure, covariate_method):
	if covariate_method == 'none':
		cluster_jxn_data_structure = add_intercept_covariate_to_cluster_jxn_data_structure(cluster_jxn_data_structure)
	return cluster_jxn_data_structure


def create_cluster_based_data_structure(tissue_specific_jxn_file, max_number_of_junctions_per_cluster, covariate_method):
	#Get raw data structure
	#Also get samples, this is the maximum possible samples after filtering
	cluster_jxn_data_structure, samples = extract_raw_cluster_jxn_data_structure(tissue_specific_jxn_file)

	# Remove clusters with more than $max_number_of_junctions_per_cluster
	cluster_jxn_data_structure = max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_number_of_junctions_per_cluster)

	# Remove samples with less than $min_reads_per_sample_in_cluster reads (summed acrosss all junctions in cluster)
	# cluster_jxn_data_structure = min_reads_per_sample_filter(cluster_jxn_data_structure, min_reads_per_sample_in_cluster)

	# Add Covariate Matrix to cluster_jxn_data_structure
	cluster_jxn_data_structure = add_covariates_to_cluster_jxn_data_structure(cluster_jxn_data_structure, covariate_method)

	return cluster_jxn_data_structure, samples


splicing_outlier_file = sys.argv[1]
tissue_specific_junction_file = sys.argv[2]
output_file = sys.argv[3]



max_number_of_junctions_per_cluster = 20
covariate_method='none'
# Extract junction data and place in compact data structure
# Keys are cluster_ids and values are jxn_counts
cluster_jxn_data_structure, all_samples = create_cluster_based_data_structure(tissue_specific_junction_file, max_number_of_junctions_per_cluster, covariate_method)

# Print to output file
f = open(splicing_outlier_file)
t = open(output_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	cluster_id = data[0]
	num_samples = len(data) - 1
	if cluster_id not in cluster_jxn_data_structure:
		print('asssumption errroro')
		pdb.set_trace()
	jxn_counts = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
	max_from_one_cluster = np.max(np.sum(jxn_counts,axis=0)/np.sum(jxn_counts))
	arr = np.asarray([max_from_one_cluster]*num_samples).astype(str)
	t.write('\t'.join(arr) + '\n')
f.close()
t.close()