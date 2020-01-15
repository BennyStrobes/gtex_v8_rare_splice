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

def get_variant_to_maf_mapping(maf_file):
	f = open(maf_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[0] + '_' + data[1] + '_' + data[4].split(':')[0] + '_' +  data[5].split(':')[0] + '_1'
		variant_id2 = data[0] + '_' + data[1] + '_' + data[5].split(':')[0] + '_' +  data[4].split(':')[0] + '_1'
		variant_id3 = data[0] + '_' + data[1] + '_' + data[4].split(':')[0] + '_' +  data[5].split(':')[0] + '_2'
		variant_id4 = data[0] + '_' + data[1] + '_' + data[5].split(':')[0] + '_' +  data[4].split(':')[0] + '_2'
		af1 = float(data[4].split(':')[1])
		af2 = float(data[5].split(':')[1])
		maf = min(af1,af2)
		# Simple error checking
		if variant_id in mapping:
			print('assumptionerror!')
			pdb.set_trace()
		mapping[variant_id] = maf
		mapping[variant_id2] = maf
		mapping[variant_id3] = maf
		mapping[variant_id4] = maf
	return mapping

outlier_file = sys.argv[1]
variant_bed_file = sys.argv[2]
maf_file = sys.argv[3]
merged_data_set_file = sys.argv[4]
merged_compressed_data_set_file = sys.argv[5]


outliers = get_outlier_calls(outlier_file)

variant_to_maf = get_variant_to_maf_mapping(maf_file)


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
	river_posterior = float(data[5])
	direction = data[6]
	pval = float(data[7])
	if test_id not in dicti:
		dicti[test_id] = ([pval], posterior, gam_posterior, river_posterior, direction)
	else:
		old_tuple = dicti[test_id]
		old_tuple[0].append(pval)
		if old_tuple[1] != posterior:
			print('assumption erro!')
			pdb.set_trace()
		if old_tuple[4] != direction:
			print('assumption erro!')
			pdb.set_trace()
		dicti[test_id] = old_tuple
f.close()
t = open(merged_compressed_data_set_file, 'w')
t.write('variant_id\tensamble_id\tamish_maf\tmedian_watershed_posterior\tmedian_gam_posterior\tmedian_river_posterior\tmedian_amish_pvalue\tmode_gtex_direction\n')
for test_id in dicti.keys():
	variant_id = test_id.split(':')[0]
	gene_id = test_id.split(':')[1]
	tupler = dicti[test_id]
	median_pval = np.median(tupler[0])

	# simple error check
	if variant_id not in variant_to_maf:
		print('assumptionerrororor')
		pdb.set_trace()
	maf = variant_to_maf[variant_id]
	if maf == 0.0:
		continue
	# Print results
	t.write(variant_id + '\t' + gene_id + '\t' + str(maf) + '\t' + str(tupler[1]) + '\t' + str(tupler[2]) + '\t' + str(tupler[3]) + '\t' + str(median_pval) + '\t' + tupler[4] + '\n')
t.close()


