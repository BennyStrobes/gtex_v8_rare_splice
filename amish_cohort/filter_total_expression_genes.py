import numpy as np 
import os
import sys
import pdb



def get_gtex_v8_genes(gtex_gene_list):
	dicti = {}
	f = open(gtex_gene_list)
	for line in f:
		line = line.rstrip()
		data = line.split()
		dicti[data[0]] = 1
	return dicti

def get_expressed_genes(count_file, tpm_file):
	f = open(count_file)
	g = open(tpm_file)
	dicti = {}
	head_count = 0
	for line in f:
		line1 = line.rstrip()
		data_count = line1.split()
		line2 = g.next().rstrip()
		data_tpm = line2.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data_tpm[0] != data_count[0]:
			print('assumption erro')
		ensamble_id = data_tpm[0]
		tpm = np.asarray(data_tpm[1:]).astype(float)
		counts = np.asarray(data_count[1:]).astype(float)
		# Compute number of samples that pass our filter
		num_samples_that_pass = len(np.where((counts >= 6) & (tpm > .1))[0])
		total_samples = len(counts)
		fraction_passed = float(num_samples_that_pass)/total_samples
		if fraction_passed >= .1:
			dicti[ensamble_id] = 1
	f.close()
	g.close()
	return dicti

count_file = sys.argv[1]
tpm_file = sys.argv[2]
log_tpm_file = sys.argv[3] # input file
gtex_gene_list = sys.argv[4]
log_tpm_filtered_file = sys.argv[5]  # output file



# Extract dictionary list of genes used in gtex v8
gtex_v8_genes = get_gtex_v8_genes(gtex_gene_list)


# Extract dictionary list of expressed genes
expressed_genes = get_expressed_genes(count_file, tpm_file)


# Stream input file and print lines that pass filter
head_count = 0
f = open(log_tpm_file)
t = open(log_tpm_filtered_file,'w')
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	ensamble_id = data[0]
	if ensamble_id in gtex_v8_genes and ensamble_id in expressed_genes:
		t.write(line + '\n')
f.close()
t.close()