import numpy as np 
import os
import sys
import pdb
import gzip


def extract_dictionary_list_of_individuals(gtex_lymphocyte_standardized_expression_file):
	valid_individuals = {}
	head_count = 0
	f = gzip.open(gtex_lymphocyte_standardized_expression_file)
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			line = line.rstrip()
			data = line.split()
			for indi_id in data[4:]:
				valid_individuals[indi_id] = 1
			continue
	f.close()
	return valid_individuals 


gtex_lymphocyte_count_file = sys.argv[1]
gtex_lymphocyte_standardized_expression_file = sys.argv[2]
filtered_gtex_lymphocyte_count_file = sys.argv[3]

# First extract dictionary list of individuals to use
# These are the individuals present in gtex_lymphocyte_standardized_expression_file and thus used in the gtex rare analysis
valid_individuals = extract_dictionary_list_of_individuals(gtex_lymphocyte_standardized_expression_file)


# Now filter gtex_lymphocyte_count_file to contain only valid_individuals
f = open(gtex_lymphocyte_count_file)
t = open(filtered_gtex_lymphocyte_count_file, 'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		indi = {}
		t.write(data[0] + '\t')
		valid_indices = []
		names = np.asarray(data[2:])
		short_names = []
		for i, ele in enumerate(names):
			indi_id = 'GTEX-' + ele.split('-')[1]
			if indi_id in valid_individuals:
				valid_indices.append(i)
				short_names.append(indi_id)
		valid_indices = np.asarray(valid_indices)
		t.write('\t'.join(short_names) + '\n')
		continue
	ensamble_id = data[0]
	counts = np.asarray(data[2:])
	filtered_counts = counts[valid_indices]
	t.write(ensamble_id + '\t' + '\t'.join(filtered_counts) + '\n')
f.close()
t.close()

