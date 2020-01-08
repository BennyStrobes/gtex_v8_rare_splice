import numpy as np 
import os
import sys
import pdb
import gzip









def get_sample_mapping(sample_mapping_file):
	dicti = {}
	data_full = np.loadtxt(sample_mapping_file, dtype=str)
	# skip header
	data = data_full[1:,:]
	num_samples = data.shape[0]
	for sample_num in range(num_samples):
		rna_seq_id = data[sample_num,0]
		vcf_id = data[sample_num,1]
		dicti[rna_seq_id] = vcf_id
	return dicti

def get_list_of_valid_individuals(individual_file):
	dicti = {}
	f = open(individual_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		dicti[data[0]] = 1
	return dicti




#######################
# Command line args
#######################
raw_count_file = sys.argv[1]  # Input file
sample_mapping_file = sys.argv[2]  # File containing sample id mapping
processed_feature_count_file = sys.argv[3]  # output file
individual_file = sys.argv[4]  # file containing individauls used in splicing outlier calling

# Dictionary list of individuals used in splicing outlier calling
valid_individuals = get_list_of_valid_individuals(individual_file)

# Create mapping from rna-seq id to vcf id
sample_mapping = get_sample_mapping(sample_mapping_file)



# Stream and print
data = np.loadtxt(raw_count_file,dtype=str,delimiter='\t')
t = open(processed_feature_count_file, 'w')

sample_names = data[0,6:]
revised_sample_names = []
valid_indices = []
for i, sample_name in enumerate(sample_names):
	rna_seq_id = sample_name.split('GM')[1].split('_')[0]
	if rna_seq_id.startswith('0'):
		rna_seq_id = rna_seq_id[1:]
	if rna_seq_id not in sample_mapping:
		print('miss')
		pdb.set_trace()
	if sample_mapping[rna_seq_id] in valid_individuals:
		revised_sample_names.append(sample_mapping[rna_seq_id])
		valid_indices.append(i)
valid_indices = np.asarray(valid_indices)
revised_sample_names = np.asarray(revised_sample_names)
head_count = 0
num_rows = data.shape[0]

for row_num in range(num_rows):
	# header
	if row_num == 0:
		t.write('\t'.join(data[0, 0:6]) + '\t' + '\t'.join(revised_sample_names) + '\n')
	else:
		filtered_counts = np.asarray(data[row_num,6:])[valid_indices]
		t.write('\t'.join(data[row_num,0:6]) + '\t' + '\t'.join(filtered_counts) + '\n')
t.close()



