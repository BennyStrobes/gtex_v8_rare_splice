import numpy as np 
import os
import sys
import pdb

# Get mapping from cluster id to ensamble id
def get_mapping_from_cluster_id_to_ensamble_id(cluster_info_file):
	mapping = {}
	head_count = 0  # To skip header
	# Stream input file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		gene_string = data[2]
		gene_arr = gene_string.split(',')
		for gene in gene_arr:
			if cluster_id not in mapping:
				mapping[cluster_id] = [gene]
			else:
				arr = mapping[cluster_id]
				arr.append(gene)
				mapping[cluster_id] = arr
	f.close()
	return mapping

def get_element_wise_minimum(pvalues1, pvalues2):
	min_pvalue = np.ones(len(pvalues1))
	for i, pvalue1 in enumerate(pvalues1):
		pvalue2 = pvalues2[i]
		if np.isnan(pvalue1):
			min_pvalue[i] = pvalue2
		elif np.isnan(pvalue2):
			min_pvalue[i] = pvalue1
		else:
			min_pvalue[i] = min(pvalue1,pvalue2)
	return min_pvalue

# Get mapping from ensamble id to min vector of pvalues
def get_mapping_from_ensamble_to_min_pvalues(cluster_level_outlier_file, cluster_to_ensambles):
	mapping = {}
	head_count = 0
	# Stream input file
	f = open(cluster_level_outlier_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			data[0] = 'GENE_ID'
			header = '\t'.join(data)
			continue
		cluster_id = data[0]
		pvalues = np.asarray(data[1:]).astype(float)
		ensamble_arr = np.unique(cluster_to_ensambles[cluster_id])
		for ensamble_id in ensamble_arr:
			# Ensmable_id never seen before
			if ensamble_id not in mapping:
				mapping[ensamble_id] = pvalues
			# Ensamble_id has been seen before
			else:
				new_pvalues = get_element_wise_minimum(pvalues, mapping[ensamble_id])
				mapping[ensamble_id] = new_pvalues
	f.close()
	return header, mapping

cluster_level_outlier_file = sys.argv[1]  # Input file
gene_level_outlier_file = sys.argv[2]  # output file
cluster_info_file = sys.argv[3]  # file containing cluster to gene mapping

# Get mapping from cluster id to ensamble id
cluster_to_ensambles = get_mapping_from_cluster_id_to_ensamble_id(cluster_info_file)


# Get mapping from ensamble id to min vector of pvalues
header, ensamble_to_pvalues = get_mapping_from_ensamble_to_min_pvalues(cluster_level_outlier_file, cluster_to_ensambles)

# Print to output file
t = open(gene_level_outlier_file, 'w')
# print header
t.write(header + '\n')
# Print line by line (gene by gene)
for ensamble_id in ensamble_to_pvalues.keys():
	t.write(ensamble_id + '\t')
	pvalues = ensamble_to_pvalues[ensamble_id]
	pvalue_arr = pvalues.astype(str)
	cleaned_pvalues = []
	for ele in pvalue_arr:
		if ele == 'nan':
			cleaned_pvalues.append('NaN')
		else:
			cleaned_pvalues.append(ele)
	t.write('\t'.join(cleaned_pvalues) + '\n')
t.close()