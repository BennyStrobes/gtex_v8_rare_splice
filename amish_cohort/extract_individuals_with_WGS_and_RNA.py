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
		dicti[ind] =1
	return dicti



vcf_file = sys.argv[1]
junction_count_input_dir = sys.argv[2]
individual_list = sys.argv[3]
individual_to_junction_file_list = sys.argv[4]



vcf_individuals = get_vcf_individuals(vcf_file)

both_individuals = {}
for file_name in os.listdir(junction_count_input_dir):
	if file_name.endswith('junc.gz') == False:
		continue
	name = file_name.split('_')[1].split('.')[0]
	name = name[2:]
	if name in vcf_individuals:
		if name not in both_individuals:
			both_individuals[name] = file_name
		else:
			print('repeated name')
			pdb.set_trace()

individuals = sorted(both_individuals.keys())
t = open(individual_list, 'w')
t2 = open(individual_to_junction_file_list, 'w')
for individual in individuals:
	t.write(individual + '\n')
	t2.write(individual + '\t' + both_individuals[individual] + '\n')
t.close()
t2.close()

