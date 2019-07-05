import numpy as np 
import os
import sys
import pdb


def get_median_stats(file_name):
	medians = []
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			# feat = data[:49]
			# pvalz = data[50:-1]
			# brain_pvalz = data[55:68]
			continue
		brain_pvalz = np.asarray(data[55:68]).astype(float)
		if sum(np.isnan(brain_pvalz) == False) >= 3:
			medians.append(np.nanmedian(brain_pvalz))
		else:
			medians.append(np.nan)
	f.close()
	return medians

def remove_bad_cols(output_file2, output_file3):
	f = open(output_file2)
	cols = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		for i, ele in enumerate(data):
			if i not in cols:
				cols[i] = []
			cols[i].append(ele)
	f.close()
	bad_indices = {}
	for col_num in cols.keys():
		if len(np.unique(cols[col_num])) < 3:
			bad_indices[col_num] =1
	f = open(output_file2)
	t = open(output_file3, 'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		for i,ele in enumerate(data):
			if i not in bad_indices:
				if i == 0:
					t.write(ele)
				else:
					t.write('\t' + ele)
		t.write('\n')
	f.close()
	t.close()

#######################
# Command line args
#######################
input_output_dir = sys.argv[1]
pvalue_threshold = sys.argv[2]

ase_file = input_output_dir + 'ase_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
splicing_file = input_output_dir + 'splicing_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
te_file = input_output_dir + 'total_expression_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'


med_ase = get_median_stats(ase_file)
med_splice = get_median_stats(splicing_file)
med_te = get_median_stats(te_file)
output_file = input_output_dir + 'median_brain_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
output_file2 = input_output_dir + 'median_observed_brain_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
output_file3 = input_output_dir + 'median_observedd_brain_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'

t = open(output_file, 'w')
t2 = open(output_file2, 'w')
f = open(ase_file)
head_count = 0
counter = 0
counter2 = 0
used = {}
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data[:49]) + '\t')
		t.write('splicing_pvalue\ttotal_expression_pvalue\tase_pvalue\tN2pair\n')
		t2.write('\t'.join(data[:49]) + '\t')
		t2.write('splicing_pvalue\ttotal_expression_pvalue\tase_pvalue\tN2pair\n')
		continue
	t.write('\t'.join(data[:49]) + '\t')
	t.write(str(med_splice[counter]) + '\t' + str(med_te[counter]) + '\t' + str(med_ase[counter]) + '\t' + data[-1] + '\n')
	if np.isnan(med_splice[counter]) == False or np.isnan(med_te[counter]) == False or np.isnan(med_ase[counter]) == False:
		if data[-1] == 'NA':
			t2.write('\t'.join(data[:49]) + '\t')
			t2.write(str(med_splice[counter]) + '\t' + str(med_te[counter]) + '\t' + str(med_ase[counter]) + '\t' + data[-1] + '\n')	
		else:
			if data[-1] not in used:
				used[data[-1]] = (data, med_splice[counter], med_te[counter], med_ase[counter])
			else:
				old_data_obj = used[data[-1]]
				old_data = old_data_obj[0]
				t2.write('\t'.join(old_data[:49]) + '\t')
				t2.write(str(old_data_obj[1]) + '\t' + str(old_data_obj[2]) + '\t' + str(old_data_obj[3]) + '\t' + str(counter2) + '\n')
				t2.write('\t'.join(data[:49]) + '\t')
				t2.write(str(med_splice[counter]) + '\t' + str(med_te[counter]) + '\t' + str(med_ase[counter]) + '\t' + str(counter2) + '\n')	
				counter2 = counter2 + 1
	counter = counter + 1
f.close()
t.close()
t2.close()

remove_bad_cols(output_file2, output_file3)
