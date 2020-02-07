import numpy as np 
import os
import sys
import pdb




def extract_data_from_a_tissue_file(tbt_file, number_of_tissues, number_of_genomic_annotations):
	head_count = 0
	genomic_annotations = []
	outlier_calls = []
	n2_pairs = []
	f = open(tbt_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields from line
		line_genomic_anno = data[0:number_of_genomic_annotations]
		line_outliers = data[number_of_genomic_annotations:(number_of_genomic_annotations+number_of_tissues)]
		line_n2_pairs = data[-1]
		# add to array
		genomic_annotations.append(np.asarray(line_genomic_anno))
		outlier_calls.append(np.asarray(line_outliers))
		n2_pairs.append(line_n2_pairs)

	f.close()
	return genomic_annotations, outlier_calls, n2_pairs

def extact_header(tbt_file, number_of_tissues, number_of_genomic_annotations):
	head_count = 0
	f = open(tbt_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			header_genomic_anno = np.asarray(data[0:number_of_genomic_annotations])
			header_outliers = np.asarray(data[number_of_genomic_annotations:(number_of_genomic_annotations+number_of_tissues)])
			header_n2_pairs = data[-1]
			continue
	f.close()
	return header_genomic_anno, header_outliers, header_n2_pairs


unsupervised_learning_input_dir = sys.argv[1]




splicing_tbt_file = sys.argv[1] + 'splicing_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
ase_tbt_file = sys.argv[1] + 'ase_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
expression_tbt_file = sys.argv[1] + 'total_expression_tbt_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'

number_of_tissues = 3
number_of_genomic_annotations = 49

splicing_genomic_annotations, splicing_outlier_calls, splicing_n2_pairs = extract_data_from_a_tissue_file(splicing_tbt_file, number_of_tissues, number_of_genomic_annotations)

ase_genomic_annotations, ase_outlier_calls, ase_n2_pairs = extract_data_from_a_tissue_file(ase_tbt_file, number_of_tissues, number_of_genomic_annotations)

expression_genomic_annotations, expression_outlier_calls, expression_n2_pairs = extract_data_from_a_tissue_file(expression_tbt_file, number_of_tissues, number_of_genomic_annotations)


# Extact header
genomic_anno_header, tissue_names, n2_pair_header = extact_header(splicing_tbt_file, number_of_tissues, number_of_genomic_annotations)

# Print to output file
output_file = unsupervised_learning_input_dir + 'tbt_cross_signal_outliers_0.01_genes_intersection_between_te_ase_splicing_features_filter_N2_pairs.txt'
t = open(output_file, 'w')
# Write header
t.write('\t'.join(genomic_anno_header))
for tissue_name in tissue_names:
	t.write('\t' + 'splicing_' + tissue_name)
for tissue_name in tissue_names:
	t.write('\t' + 'expression_' + tissue_name)
for tissue_name in tissue_names:
	t.write('\t' + 'ase_' + tissue_name)
t.write('\t' + n2_pair_header + '\n')
# Write lines to output file
num_lines = len(splicing_n2_pairs)
for line_number in range(num_lines):
	t.write('\t'.join(splicing_genomic_annotations[line_number]) + '\t' + '\t'.join(splicing_outlier_calls[line_number]) + '\t' + '\t'.join(expression_outlier_calls[line_number]) + '\t' + '\t'.join(ase_outlier_calls[line_number]) + '\t' + expression_n2_pairs[line_number] + '\n')
t.close()

# Create connection file
output_file = unsupervised_learning_input_dir + 'tbt_cross_signal_outliers_0.01_genes_intersection_between_te_ase_splicing_connection_file.txt'
total_nodes = number_of_tissues*number_of_tissues
t = open(output_file,'w')
for ii in range(total_nodes):
	for jj in range(total_nodes):
		if jj != 0:
			t.write('\t')
		if ii == jj:
			t.write('0')
		else:
			# if in same outlier signal, place edge
			if ii/number_of_tissues == jj/number_of_tissues:
				t.write('1')
			elif ii%number_of_tissues == jj%number_of_tissues: # If from same tissue, place edge
				t.write('1')
			else:
				t.write('0')
	t.write('\n')
t.close()

