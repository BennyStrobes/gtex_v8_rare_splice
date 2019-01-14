import numpy as np 
import os
import sys
import pdb
import pickle



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
		tissue_splicing_outlier_file = splicing_outlier_dir + tissue + '_covariate_method_' + covariate_method + '_merged_emperical_pvalue_gene_level.txt'
		# Fill in multi_tissue_outliers for this tissue
		multi_tissue_outliers = update_multi_tissue_outliers_for_one_tissue(tissue_splicing_outlier_file, multi_tissue_outliers)
	return multi_tissue_outliers

# Create object containing mapping from individual to number of global_outliers
def count_number_of_multi_tissue_outliers_per_individual(multi_tissue_outliers, min_number_of_expressed_tissues, pvalue_threshold):
	# Initialize object
	individual_to_number_outliers = {}
	# Loop through clusters
	for cluster_id in multi_tissue_outliers.keys():
		# Loop through individuals for this cluster
		for individual in multi_tissue_outliers[cluster_id].keys():
			# Add individual to object if not already there 
			if individual not in individual_to_number_outliers:
				individual_to_number_outliers[individual] = 0
			# Vector of pvalues across tissues
			multi_tissue_pvalue_vector = multi_tissue_outliers[cluster_id][individual]
			# Ignore (individual, cluster pairs) that are expressed in fewer than $min_number_of_expressed_tissues
			if len(multi_tissue_pvalue_vector) < min_number_of_expressed_tissues:
				continue
			# Check if individual is a multi tissue outlier for this cluster
			if np.median(multi_tissue_pvalue_vector) < pvalue_threshold:
				individual_to_number_outliers[individual] = individual_to_number_outliers[individual] + 1
	return individual_to_number_outliers

# Print mapping of individual id to number of global outliers to output file
def print_individual_to_number_outliers(individual_to_number_outliers, output_file):
	# Create array of tuples from dictionary
	outliers = []
	for individual in individual_to_number_outliers.keys():
		num_outliers = individual_to_number_outliers[individual]
		outliers.append((individual, num_outliers))

	# Sort individuals by number of multitissue outliers
	outliers.sort(key=lambda x: x[1],reverse=True)

	# Open output file
	t = open(output_file, 'w')
	# Print header of output file
	t.write('individual\tnumber_of_multi_tissue_outliers\n')
	for outlier_tuple in outliers:
		t.write(outlier_tuple[0] + '\t' + str(outlier_tuple[1]) + '\n')
	t.close()

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

# Print mapping of individual id to number of global outliers to output file while limiting to only european ancestry individuals
def print_individual_to_number_outliers_european_only(individual_to_number_outliers_output_file, individual_to_number_outliers_ea_only_output_file, ea_individuals):
	# To identify header
	head_count = 0
	# Open input file handle
	f = open(individual_to_number_outliers_output_file)
	# open output file handle
	t = open(individual_to_number_outliers_ea_only_output_file, 'w')
	# Stream input file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		# standard line
		# Skip non european ancestry individuals
		if data[0] not in ea_individuals:
			continue
		t.write(line + '\n')
	f.close()
	t.close()

# Extract dictionary list of individuals that are "good" (to keep) samples
def extract_non_global_outlier_individuals(individual_to_number_outliers_output_file):
	counts = np.loadtxt(individual_to_number_outliers_output_file,dtype=str,delimiter='\t')[1:,1].astype(float)
	max_number_global_outliers = np.percentile(counts,75) + 3.0*(np.percentile(counts,75) - np.percentile(counts,25))
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
		if num_outliers <= max_number_global_outliers:
			individuals_to_keep[individual] = 1
	f.close()
	return individuals_to_keep

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
		# And print line
		t.write('\t'.join(filtered_data) + '\n')
	f.close()
	t.close()

# For each tissue specific outlier file, remove individuals that are global outliers
def remove_global_outliers_from_tissue_specific_outlier_calls(individual_to_number_outliers_output_file, tissue_names, covariate_method, splicing_outlier_dir, suffix):
	# Extract dictionary list of individuals that are "good" (to keep) samples
	individuals_to_keep = extract_non_global_outlier_individuals(individual_to_number_outliers_output_file)
	# Loop through tissues
	for tissue_name in tissue_names:
		# Splicing outlier input file
		outlier_input_file = splicing_outlier_dir + tissue_name + '_covariate_method_' + covariate_method + '_merged_' + suffix + '.txt'
		# Filtered splicing outlier output file
		outlier_output_file = splicing_outlier_dir + tissue_name + '_covariate_method_' + covariate_method + '_no_global_outliers_merged_' + suffix + '.txt'
		# Filter individuals (columns) in outlier file
		filter_individuals_of_outlier_file(outlier_input_file, outlier_output_file, individuals_to_keep)

# For each tissue specific outlier file, remove individuals that are not european ancestry
def remove_non_european_ancestry_individuals_from_tissue_specific_outlier_calls(ea_individuals, tissue_names, covariate_method, splicing_outlier_dir, suffix):
	# Loop through tissues
	for tissue_name in tissue_names:
		# Splicing outlier input file
		outlier_input_file = splicing_outlier_dir + tissue_name + '_covariate_method_' + covariate_method + '_no_global_outliers_merged_' + suffix + '.txt'
		# Filtered splicing outlier output file
		outlier_output_file = splicing_outlier_dir + tissue_name + '_covariate_method_' + covariate_method + '_no_global_outliers_ea_only_merged_' + suffix + '.txt'
		# Filter individuals (columns) in outlier file
		filter_individuals_of_outlier_file(outlier_input_file, outlier_output_file, ea_individuals)

# File containing list of gtex tissue names
tissue_names_file = sys.argv[1]
# Directory containing splicing outlier results
splicing_outlier_dir = sys.argv[2]
# How covariates were handled in outlier calling (used for filenames in this script)
covariate_method = sys.argv[3]
# To be eligible to be a "global outlier", (individual, cluster) pair must be expressed in at least $min_number_of_expressed_tissues
min_number_of_expressed_tissues = int(sys.argv[4])
# Anything less than p=$pvalue_threshold is considered an outlier
pvalue_threshold = float(sys.argv[5])
# File containing list of individuals of european ancestry
european_ancestry_individual_list = sys.argv[6]


# Extract vector of tissue names
tissue_names = extract_tissue_names(tissue_names_file)

# Create object that maps from [gene_id][donor_id] to a vector of pvalues (each element corresponding to one tissue)
multi_tissue_outliers = create_multi_tissue_outlier_object(tissue_names, splicing_outlier_dir, covariate_method)

# Create object containing mapping from individual to number of global_outliers
individual_to_number_outliers = count_number_of_multi_tissue_outliers_per_individual(multi_tissue_outliers, min_number_of_expressed_tissues, pvalue_threshold)

# Print mapping of individual id to number of global outliers to output file
individual_to_number_outliers_output_file = splicing_outlier_dir + 'number_of_multi_tissue_outliers.txt'
print_individual_to_number_outliers(individual_to_number_outliers, individual_to_number_outliers_output_file)

# Extract dictionary list of european ancestry individuals
ea_individuals = extract_dictionary_list_of_european_ancestry_individuals(european_ancestry_individual_list)

# Print mapping of individual id to number of global outliers to output file while limiting to only european ancestry individuals
individual_to_number_outliers_ea_only_output_file = splicing_outlier_dir + 'number_of_multi_tissue_outliers_ea_only.txt'
print_individual_to_number_outliers_european_only(individual_to_number_outliers_output_file, individual_to_number_outliers_ea_only_output_file, ea_individuals)

# For each tissue specific outlier file, remove individuals that are global outliers (for cluster level outlier file)
remove_global_outliers_from_tissue_specific_outlier_calls(individual_to_number_outliers_output_file, tissue_names, covariate_method, splicing_outlier_dir, 'emperical_pvalue')

# For each tissue specific outlier file, remove individuals that are not european ancestry (for cluster level outlier file)
remove_non_european_ancestry_individuals_from_tissue_specific_outlier_calls(ea_individuals, tissue_names, covariate_method, splicing_outlier_dir, 'emperical_pvalue')

# For each tissue specific outlier file, remove individuals that are global outliers (for gene level outlier file)
remove_global_outliers_from_tissue_specific_outlier_calls(individual_to_number_outliers_output_file, tissue_names, covariate_method, splicing_outlier_dir, 'emperical_pvalue_gene_level')

# For each tissue specific outlier file, remove individuals that are not european ancestry (for gene level outlier file)
remove_non_european_ancestry_individuals_from_tissue_specific_outlier_calls(ea_individuals, tissue_names, covariate_method, splicing_outlier_dir, 'emperical_pvalue_gene_level')