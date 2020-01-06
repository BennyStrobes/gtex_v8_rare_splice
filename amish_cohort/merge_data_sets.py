import numpy as np 
import os
import sys
import pdb




def get_outlier_calls(splicing_outlier_file):
	f = open(splicing_outlier_file)
	head_count = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			samples = np.asarray(data[1:])
			continue
		ensamble_id = data[0].split('.')[0]
		pvalz = np.asarray(data[1:])
		for i, pval in enumerate(pvalz):
			sample_name = samples[i]
			if pval == 'NA':
				continue
			dicti[sample_name + '_' + ensamble_id] = pval
	f.close()
	return dicti


outlier_file = sys.argv[1]
variant_bed_file = sys.argv[2]
merged_data_set_file = sys.argv[3]
merged_compressed_data_set_file = sys.argv[4]


outliers = get_outlier_calls(outlier_file)


f = open(variant_bed_file)
t = open(merged_data_set_file, 'w')
head_count = 0

for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\tamish_pvalue\n')
		continue
	sample_name = data[2]
	ensamble_id = data[1].split('.')[0]
	test_name = sample_name + '_' + ensamble_id
	if test_name not in outliers:
		continue
	t.write(line + '\t' + outliers[test_name] + '\n')
f.close()
t.close()

dicti = {}
head_count = 0
f = open(merged_data_set_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	test_id = data[0] + ':' + data[1]
	posterior = float(data[3])
	gam_posterior = float(data[4])
	pval = float(data[5])
	if test_id not in dicti:
		dicti[test_id] = ([pval], posterior, gam_posterior)
	else:
		old_tuple = dicti[test_id]
		old_tuple[0].append(pval)
		if old_tuple[1] != posterior:
			print('assumption erro!')
			pdb.set_trace()
		dicti[test_id] = old_tuple
f.close()
t = open(merged_compressed_data_set_file, 'w')
t.write('variant_id\tensamble_id\tmedian_watershed_posterior\tmedian_gam_posterior\tmedian_amish_pvalue\n')
for test_id in dicti.keys():
	variant_id = test_id.split(':')[0]
	gene_id = test_id.split(':')[1]
	tupler = dicti[test_id]
	median_pval = np.median(tupler[0])
	t.write(variant_id + '\t' + gene_id + '\t' + str(tupler[1]) + '\t' + str(tupler[2]) + '\t' + str(median_pval) + '\n')
t.close()


