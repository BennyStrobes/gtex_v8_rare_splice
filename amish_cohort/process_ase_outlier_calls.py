import numpy as np 
import os
import sys
import pdb



def get_sample_mapping(sample_mapping_file):
	dicti = {}
	data_full = np.loadtxt(sample_mapping_file, dtype=str)
	# skip header
	data = data_full[1:,:]
	num_samples = data.shape[0]
	for sample_num in range(num_samples):
		rna_seq_id = data[sample_num,0]
		vcf_id = data[sample_num,1]
		dicti['GM' + rna_seq_id] = vcf_id
	return dicti







ase_outlier_calls_file = sys.argv[1]
processed_ase_outlier_calls_file = sys.argv[2]
sample_mapping_file = sys.argv[3]


sample_mapping = get_sample_mapping(sample_mapping_file)

f = open(ase_outlier_calls_file)
t = open(processed_ase_outlier_calls_file, 'w')


head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(data[0])
		for ele in data[1:]:
			if ele.startswith('GM0'):
				word = 'GM' + ele[3:]
			else:
				word = ele
			if word not in sample_mapping:
				pdb.set_trace()
			else:
				t.write('\t' + sample_mapping[word])
		t.write('\n')
		continue
	t.write(line + '\n')
f.close()
t.close()