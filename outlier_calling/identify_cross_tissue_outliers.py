import numpy as np 
import os
import sys
import pdb



# Extract vector of tissue names
def extract_tissue_names(tissue_names_file):
	# Initialze output vector
	tissue_names = []
	# Stream file
	f = open(tissue_names_file)
	for line in f:
		line = line.rstrip()
		tissue_names.append(line)
	f.close()
	return np.asarray(tissue_names)

# Extract dictionary list of european ancestry individuals
def extract_dictionary_list_of_european_ancestry_individuals(european_ancestry_individual_list):
	# Initialize dictionary list
	ea_individuals = {}
	# Stream input file
	f = open(european_ancestry_individual_list)
	for line in f:
		line = line.rstrip()
		data = line.split()
		ea_individuals[line] = 1
	f.close()
	return ea_individuals

def update_multi_tissue_outliers_for_one_tissue(tissue_splicing_outlier_file, multi_tissue_outliers):
	# Variable used to identify header
	head_count = 0
	# Stream outlier file
	f = open(tissue_splicing_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# HEADER
		if head_count == 0:
			head_count = head_count + 1
			# Ordered list of individuals for this tissues
			individuals = np.asarray(data[1:])
			continue
		# Parse standard line
		cluster_id = data[0]
		pvalues = np.asarray(data[1:]).astype(float)
		# If cluster_id has never been seen before, add to multi_tissue_outliers object
		if cluster_id not in multi_tissue_outliers:
			multi_tissue_outliers[cluster_id] = {}
		# Loop through individuals
		for i, individual in enumerate(individuals):
			# pvalue for this individual
			pvalue = pvalues[i]
			# Now add this individual to object
			if individual not in multi_tissue_outliers[cluster_id]:
				multi_tissue_outliers[cluster_id][individual] = []
			multi_tissue_outliers[cluster_id][individual].append(pvalue)
	f.close()
	return multi_tissue_outliers


# Create object that maps from [cluster_id][donor_id] to a vector of pvalues (each element corresponding to one tissue)
def create_multi_tissue_outlier_object(tissue_names, splicing_outlier_dir, covariate_method):
	# Initialize object
	multi_tissue_outliers = {}
	# Fill in multi_tissue_outliers for each tissue in series
	for tissue in tissue_names:
		# outlier file for this tissue
		tissue_splicing_outlier_file = splicing_outlier_dir + tissue + '_covariate_method_' + covariate_method + '_merged_emperical_pvalue.txt'
		# Fill in multi_tissue_outliers for this tissue
		multi_tissue_outliers = update_multi_tissue_outliers_for_one_tissue(tissue_splicing_outlier_file, multi_tissue_outliers)
	return multi_tissue_outliers

def extract_all_individuals(tissue_names, splicing_outlier_dir, covariate_method):
	individuals = []
	for tissue in tissue_names:
		tissue_splicing_outlier_file = splicing_outlier_dir + tissue + '_covariate_method_' + covariate_method + '_merged_emperical_pvalue.txt'
		f = open(tissue_splicing_outlier_file)
		head_count = 0
		for line in f:
			if head_count == 0:
				head_count = head_count + 1
				line = line.rstrip()
				data = line.split()
				for ele in data[1:]:
					individuals.append(ele)
				continue
		f.close()
	return sorted(np.unique(individuals))

# Print median(pvalue) outlier file
def print_cross_tissue_outlier_file(all_individuals, multi_tissue_outliers, min_number_of_expressed_tissues, cross_tissue_outlier_file):
	# Open output file handle
	t = open(cross_tissue_outlier_file, 'w')
	# Print header
	t.write('CLUSTER_ID\t' + '\t'.join(all_individuals) + '\n')
	clusters = multi_tissue_outliers.keys()
	# Loop through clusters (row)
	for cluster_id in clusters:
		t.write(cluster_id)
		# Loop through individuals (columns)
		for individual in all_individuals:
			# Individual never seen for this cluster
			if individual not in multi_tissue_outliers[cluster_id]:
				t.write('\tNaN')
				continue
			# Too few tissues for (individual, cluster) pair
			if len(multi_tissue_outliers[cluster_id][individual]) < min_number_of_expressed_tissues:
				t.write('\tNaN')
				continue
			# Extract median(pvalue)
			median_pvalue = np.median(multi_tissue_outliers[cluster_id][individual])
			# Print median(pvalue)
			t.write('\t' + str(median_pvalue))
		t.write('\n')
	t.close()

def line_is_all_nan(filtered_data):
	binary = True
	for ele in filtered_data[1:]:
		if ele != 'NaN':
			binary = False
	return binary

# Filter individuals (columns) in outlier file
def filter_individuals_of_outlier_file(outlier_input_file, outlier_output_file, individuals_to_keep):
	# For header 
	head_count =  0
	# open output file handle
	t = open(outlier_output_file, 'w')
	# Stream input file handle
	f = open(outlier_input_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			# Learn column indices to be kept
			good_columns = []
			for i, individual_id in enumerate(data):
				if i == 0 or individual_id in individuals_to_keep:
					good_columns.append(i)
			# Filter data to only those columns
			filtered_data = np.asarray(data)[good_columns]
			# And print header
			t.write('\t'.join(filtered_data) + '\n')
			continue
		# Filter data to only outlier scores from ea individuals
		filtered_data = np.asarray(data)[good_columns]
		# Check to make sure line isn't all NaNs
		if line_is_all_nan(filtered_data):
			continue
		# And print line
		t.write('\t'.join(filtered_data) + '\n')
	f.close()
	t.close()


# Extract dictionary list of individuals that are "good" (to keep) samples
def extract_non_global_outlier_individuals(individual_to_number_outliers_output_file, num_global_outlier_clusters):
	# Initialize list
	individuals_to_keep = {}
	# For header
	head_count = 0
	# Stream input file
	f = open(individual_to_number_outliers_output_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Standard line
		num_outliers = int(data[1])
		individual = data[0]
		if num_outliers < num_global_outlier_clusters:
			individuals_to_keep[individual] = 1
	f.close()
	return individuals_to_keep


# File containing list of gtex tissue names
tissue_names_file = sys.argv[1]
# Directory containing splicing outlier results
splicing_outlier_dir = sys.argv[2]
# How covariates were handled in outlier calling (used for filenames in this script)
covariate_method = sys.argv[3]
# To be eligible to be a "global outlier", (individual, cluster) pair must be expressed in at least $min_number_of_expressed_tissues
min_number_of_expressed_tissues = int(sys.argv[4])
# File containing list of individuals of european ancestry
european_ancestry_individual_list = sys.argv[5]
# A donor is called a "global outlier" if it is a median(pvalue) outlier in at least $num_global_outlier_clusters clusters
num_global_outlier_clusters = int(sys.argv[6])


# Extract vector of tissue names
tissue_names = extract_tissue_names(tissue_names_file)

# Vector of all individuals (including non-EA ancestry and global outliers)
all_individuals = extract_all_individuals(tissue_names, splicing_outlier_dir, covariate_method)

# Extract dictionary list of european ancestry individuals
ea_individuals = extract_dictionary_list_of_european_ancestry_individuals(european_ancestry_individual_list)


# Create object that maps from [cluster_id][donor_id] to a vector of pvalues (each element corresponding to one tissue)
multi_tissue_outliers = create_multi_tissue_outlier_object(tissue_names, splicing_outlier_dir, covariate_method)

# Print median(pvalue) outlier file
cross_tissue_outlier_file = splicing_outlier_dir + 'cross_tissue_covariate_method_none_emperical_pvalue.txt'
print_cross_tissue_outlier_file(all_individuals, multi_tissue_outliers, min_number_of_expressed_tissues, cross_tissue_outlier_file)

# Extract dictionary list of individuals that are "good" (to keep; non-global outliers) samples
individual_to_number_outliers_output_file = splicing_outlier_dir + 'number_of_multi_tissue_outliers.txt'
individuals_non_global_outliers = extract_non_global_outlier_individuals(individual_to_number_outliers_output_file, num_global_outlier_clusters)

# Filter individuals (columns) in outlier file
cross_tissue_outlier_no_global_file = splicing_outlier_dir + 'cross_tissue_covariate_method_none_no_global_emperical_pvalue.txt'
filter_individuals_of_outlier_file(cross_tissue_outlier_file, cross_tissue_outlier_no_global_file, individuals_non_global_outliers)

cross_tissue_outlier_no_global_no_ea_file = splicing_outlier_dir + 'cross_tissue_covariate_method_none_no_global_ea_only_emperical_pvalue.txt'
filter_individuals_of_outlier_file(cross_tissue_outlier_no_global_file, cross_tissue_outlier_no_global_no_ea_file, ea_individuals)


