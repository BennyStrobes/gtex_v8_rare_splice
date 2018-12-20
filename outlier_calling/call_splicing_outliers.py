import numpy as np 
import os
import sys
import pdb
import pickle
import time
import pystan
import dirichlet_multinomial_glm as dm_glm





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

# For parallelization purposes
def parallelization_start_and_end(num_tasks, job_number, total_jobs):
	tasks_per_job = (num_tasks/total_jobs) + 1
	start_task = job_number*tasks_per_job
	end_task = (job_number + 1)*tasks_per_job -1 
	return start_task, end_task


# Print outlier calling dm results to output file
def outlier_calling_print_helper(arr, cluster_samples, all_samples, t, cluster_id):
	counter = 0
	# Print Row id
	t.write(cluster_id)
	# Loop through all samples
	for sample in all_samples:
		# If sample in cluster specific samples
		if sample in cluster_samples:
			ele = str(arr[counter])
			counter = counter + 1
		# Samples was filtered out of this cluster
		else:
			ele = 'NaN'
		t.write('\t' + ele)
	t.write('\n')
	if counter != len(cluster_samples):
		print('print helper error!')
		pdb.set_trace()
	t.flush()
	return t



# Call splicing outliers for each cluster
# Also write to output
def call_splicing_outliers_shell(output_root, cluster_jxn_data_structure, all_samples, start_number, end_number):
	# Pystan optimizizer
	#DM_GLM = pystan.StanModel(file = "dm_glm_multi_conc.stan")
	DM_GLM = pickle.load(open('/home-1/bstrobe1@jhu.edu/scratch/gtex_v8/rare_var/gtex_v8_rare_splice/outlier_calling/dm_glm_multi_conc.pkl', 'rb'))

	#Initialize output files
	t_MD = open(output_root + '_md.txt','w')  # Filehandle for matrix of mahalanobis distances
	t_pvalue = open(output_root + '_emperical_pvalue.txt', 'w')  # Filehandle for matrix of pvalues
	# Write headers for output files
	t_MD.write('CLUSTER_ID\t' + '\t'.join(all_samples) + '\n')
	t_pvalue.write('CLUSTER_ID\t' + '\t'.join(all_samples) + '\n')

	start_time = time.time()
	county = 0 
	total = 0
	# Loop through cluster_id
	for counter, cluster_id in enumerate(sorted(cluster_jxn_data_structure.keys())):
		# Skip cluster_ids not in this parallelization run
		if counter < start_number or counter > end_number:
			continue
		####################################################################
		# Timing/Debugging related stuff
		####################################################################
		curr_time = time.time()
		diff = curr_time-start_time
		start_time = curr_time
		time_min = diff/60.0
		print(time_min)
		print('##################')
		print('START ' + cluster_id)


		####################################################################
		# Actual Analysis
		####################################################################
		# Extract jxn matrix for this gene
		X = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# Extract sample ids used in THIS cluster (NOTE: Different from all_samples)
		cluster_samples = cluster_jxn_data_structure[cluster_id]['samples']
		# Extract covariate matrix (dim len(cluster_samples)Xnum_cov)
		cov_mat = cluster_jxn_data_structure[cluster_id]['covariate_matrix']
		# Run outlier analysis:
		# Return:
		#   1: mahalanobis_distances: vector length num_samples where each element is the mahalanobis distance for that sample
		#   2. pvalues: vector of length num_samples where each element is the pvalue for that sample
		#   3. alpha: vector of length num_jxns which defines the fitted dirichlet multinomial distribution
		try:
			mahalanobis_distances, pvalues, alpha = dm_glm.run_dm_outlier_analysis(X, cov_mat, DM_GLM)


			####################################################################
			# Print results to output file
			####################################################################
			# Print Mahalanobis distance results to output file
			t_MD = outlier_calling_print_helper(mahalanobis_distances, cluster_samples, all_samples, t_MD, cluster_id)
			# Print emperical pvalue resutls to output file
			t_pvalue = outlier_calling_print_helper(pvalues, cluster_samples, all_samples, t_pvalue, cluster_id)
		except:
			print('miss: ' + cluster_id)
	print(county)
	print(total)
	t_MD.close()
	t_pvalue.close()




########################
# Command Line args
########################
# Name of GTEx tissue
tissue_type = sys.argv[1]
# Filename containing junction counts for current tissue
tissue_specific_jxn_file = sys.argv[2]
# What type of covariates to include
covariate_method = sys.argv[3]
# Throw out clusters with more than $max_number_of_junctions_per_cluster
max_number_of_junctions_per_cluster = int(sys.argv[4])
# Stem used for output files
output_root = sys.argv[5]
# Following are used for parallelization purposes
job_number = int(sys.argv[6])
total_jobs = int(sys.argv[7])


# Extract junction data and place in compact data structure
# Keys are cluster_ids and values are jxn_counts
cluster_jxn_data_structure, all_samples = create_cluster_based_data_structure(tissue_specific_jxn_file, max_number_of_junctions_per_cluster, covariate_method)


########################################
# To skip running above line
#f = open(output_root + 'cluster_data_struct.pkl', 'wb')
#pickle.dump(cluster_jxn_data_structure, f)
#f.close()

#f = open(output_root + 'all_samples.pkl', 'wb')
#pickle.dump(all_samples, f)
#f.close()

#cluster_jxn_data_structure = pickle.load(open(output_root + 'cluster_data_struct.pkl','rb'))
#all_samples = pickle.load(open(output_root + 'all_samples.pkl','rb'))
########################################

#For parallelization purposes
start_number, end_number = parallelization_start_and_end(len(cluster_jxn_data_structure), job_number, total_jobs)

# Call splicing outliers for each cluster
# Also write to output
call_splicing_outliers_shell(output_root, cluster_jxn_data_structure, all_samples, start_number, end_number)
