import numpy as np 
import os
import sys
import pdb
import gzip


def get_vcf_individuals(vcf_file):
	f = gzip.open(vcf_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('#') == False:
			break
		if line.startswith('#CHROM'):
			data = line.split('\t')
			indi = data[9:]
	dicti = {}
	for ind in indi:
		dicti[ind] = 0
	return dicti


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


vcf_file = sys.argv[1]
junction_count_input_dir = sys.argv[2]
sample_mapping_file = sys.argv[3]
individual_list = sys.argv[4]
individual_to_junction_file_list = sys.argv[5]



vcf_individuals = get_vcf_individuals(vcf_file)

sample_mapping = get_sample_mapping(sample_mapping_file)
used = {}
both_individuals = {}
for file_name in os.listdir(junction_count_input_dir):
	if file_name.endswith('junc.gz') == False:
		continue
	rna_seq_id = file_name.split('_')[0].split('GM')[1]
	if rna_seq_id[0] == '0':
		pdb.set_trace()
		rna_seq_id = rna_seq_id[1:]
	if rna_seq_id in used:
		pdb.set_trace()
	if rna_seq_id in sample_mapping:
		vcf_id = sample_mapping[rna_seq_id]
		if vcf_id in vcf_individuals:
			if vcf_id not in both_individuals:
				vcf_individuals[vcf_id] = 1
				both_individuals[vcf_id] = junction_count_input_dir + file_name
				used[rna_seq_id] = 1
			else:
				print('repeated name')
				pdb.set_trace()
		else:
			print('miss ' + rna_seq_id)
	else:
		pdb.set_trace()

individuals = sorted(both_individuals.keys())
t = open(individual_list, 'w')
t2 = open(individual_to_junction_file_list, 'w')
for individual in individuals:
	t.write(individual + '\n')
	t2.write(individual + '\t' + both_individuals[individual] + '\n')
t.close()
t2.close()


