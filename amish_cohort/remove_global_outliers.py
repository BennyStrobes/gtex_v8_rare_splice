import numpy as np 
import os
import sys
import pdb



def create_mapping_from_individual_to_num_outliers(input_file, threshold):
	mapping = {}
	head_count = 0
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			# Array of individual ids
			indi_arr = np.asarray(data[1:])
			for indi in indi_arr:
				mapping[indi] = 0
			continue
		gene_id = data[0]
		pvalz = np.asarray(data[1:]).astype(float)
		for i, pval in enumerate(pvalz):
			indi = indi_arr[i]
			if pval < threshold:
				mapping[indi] = mapping[indi] + 1
	f.close()
	num_outlier_vec = []
	for indi in mapping.keys():
		num_outlier_vec.append(mapping[indi])
	return mapping, np.asarray(num_outlier_vec)

def extract_non_global_outlier_individuals(individual_to_num_outliers, counts):
	non_global_outliers = {}
	max_number_global_outliers = np.percentile(counts,75) + 1.5*(np.percentile(counts,75) - np.percentile(counts,25))
	for indi in individual_to_num_outliers.keys():
		if individual_to_num_outliers[indi] <= max_number_global_outliers:
			non_global_outliers[indi] = 1
	return non_global_outliers

def print_outlier_file_without_global_outliers(input_file, output_file, individuals_to_keep):
	# For header 
	head_count =  0
	# open output file handle
	t = open(output_file, 'w')
	# Stream input file handle
	f = open(input_file)
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

#####################
# Command line args
#####################
# Outlier file
input_file = sys.argv[1]
# Outlier file with no global outliers
output_file = sys.argv[2]


# Create mapping from individual id to number of outliers
individual_to_num_outliers, num_outlier_vec = create_mapping_from_individual_to_num_outliers(input_file, .0027)

# Extract non-global-outlier individuals
non_global_outliers = extract_non_global_outlier_individuals(individual_to_num_outliers, num_outlier_vec)

# Print outlier file filtering out global outliers
print_outlier_file_without_global_outliers(input_file, output_file, non_global_outliers)