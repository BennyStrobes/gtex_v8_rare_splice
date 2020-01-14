import numpy as np 
import os
import sys
import pdb


def extract_gtex_count_data(filtered_gtex_lymphocyte_count_file):
	f = open(filtered_gtex_lymphocyte_count_file)
	head_count = 0
	ensamble_to_gtex_counts = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			gtex_sample_ids = np.asarray(data[1:])
			continue
		ensamble_id = data[0]
		counts = np.asarray(data[1:])
		# Simple error checking
		if ensamble_id in ensamble_to_gtex_counts:
			print('assumption error!')
			pdb.set_trace()
		ensamble_to_gtex_counts[ensamble_id] = counts
	f.close()
	return ensamble_to_gtex_counts, gtex_sample_ids

def extract_gtex_gene_length_data(gtex_gene_length_file):
	f = open(gtex_gene_length_file)
	head_count = 0
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[0]
		gene_length = data[1]
		mapping[ensamble_id] = gene_length
	f.close()
	return mapping

# Gtex count file
filtered_gtex_lymphocyte_count_file = sys.argv[1]
# amish count file
processed_feature_count_file = sys.argv[2]
# Merged (concatenated count file)
concatenated_processed_feature_count_file = sys.argv[3]
# Gene length file for gtex
gtex_gene_length_file = sys.argv[4]



# Create mapping from ensamble id to gtex gene length
ensamble_to_gtex_gene_length = extract_gtex_gene_length_data(gtex_gene_length_file)

# Create mapping from ensamble_id to gtex counts (across samples)
ensamble_to_gtex_counts, gtex_sample_names = extract_gtex_count_data(filtered_gtex_lymphocyte_count_file)

# MERGE TWO FILES
f = open(processed_feature_count_file)
t = open(concatenated_processed_feature_count_file, 'w')

head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data[0:6]) + '\tgtex_Length\t' + '\t'.join(data[6:]) + '\t' + '\t'.join(gtex_sample_names) + '\n')
		continue
	ensamble_id = data[0]
	if ensamble_id not in ensamble_to_gtex_counts:
		continue
	if ensamble_id not in ensamble_to_gtex_gene_length:
		print('miss')
		continue
	gtex_counts = ensamble_to_gtex_counts[ensamble_id]
	gtex_length = ensamble_to_gtex_gene_length[ensamble_id]
	t.write('\t'.join(data[0:6]) + '\t' + gtex_length + '\t' + '\t'.join(data[6:]) + '\t' + '\t'.join(gtex_counts) + '\n')
f.close()
t.close()