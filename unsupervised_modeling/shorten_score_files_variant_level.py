import numpy as np 
import os
import sys
import pdb
import random
import gzip



def get_categorical(valler):
	prob = -np.log10(float(valler) + .000001)
	if prob < 1:
		val = 1
	elif prob >= 1 and prob < 4:
		val = 2
	elif prob >= 4:
		val = 3
	else:
		print('assumption eroror')
		pdb.set_trace()
	return val











score_dir = sys.argv[1]


full_new_file = score_dir + 'fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt.gz'
full_old_file = score_dir + 'old_fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors.txt'

shortened_new_file = score_dir + 'fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors_short.txt'
shortened_old_file = score_dir + 'old_fully_observed_te_ase_splicing_outliers_gene_pvalue_0.01_outlier_fraction_.01_pseudocount_30_exact_inference_apply_to_all_variants_posteriors_short.txt'

output_file = score_dir + 'variant_level_comparison.txt'
# Randomly shorten new file and keep track of tests used
if full_new_file.endswith('.gz'):
	f = gzip.open(full_new_file)
else:
	f = open(full_new_file)
t = open(shortened_new_file,'w')
tests = {}
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	if random.uniform(0,1) < .005:
		t.write(line + '\n')
		tests[data[0]] = 1
f.close()
t.close()
# Shorten old file to be the same tests
f = open(full_old_file)
t = open(shortened_old_file, 'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	if data[0] in tests:
		t.write(line + '\n')
f.close()
t.close()

# Merge files
outliers = {}
scores = {}
f = open(shortened_new_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	sample_name = data[0]
	outliers[sample_name] = data[1:4]
	scores[sample_name] = data[7:]
f.close()
head_count = 0
f = open(shortened_old_file)
t = open(output_file, 'w')
t.write('ase_change\tsplicing_outlier\tte_outlier\tase_outlier\tsplicing_river_old\tte_river_old\tase_river_old\tsplicing_watershed_old\tte_watershed_old\tase_watershed_old\t')
t.write('splicing_outlier_new\tte_outlier_new\tase_outlier_new\tsplicing_river\tte_river\tase_river\tsplicing_watershed\tte_watershed\tase_watershed\n')
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	sample_name = data[0]
	if sample_name not in outliers:
		continue
	line_outlier = data[1:4]
	line_score = data[7:]
	switched = 'no_change'
	if line_outlier[2] == 'NaN' or outliers[sample_name][2] == 'NaN':
		if line_outlier[2] == outliers[sample_name][2]:
			switched = 'no_change'
		elif line_outlier[2] == 'NaN' and outliers[sample_name][2] != 'NaN':
			switched = 'NaN->observed'
		elif line_outlier[2] != 'NaN' and outliers[sample_name][2] == 'NaN':
			switched = 'observed->NaN'
		else:
			print('assumption eroro')
			pdb.set_trace()
	else:
		cat_p_1 = get_categorical(line_outlier[2])
		cat_p_2 = get_categorical(outliers[sample_name][2])
		if cat_p_1 != cat_p_2:
			if cat_p_1 > cat_p_2:
				switched = 'weaker_outlier'
			elif cat_p_1 < cat_p_2:
				switched = 'stronger_outlier'
			else:
				print('assumptine eroro')
				pdb.set_trace()
	if line_outlier[0] != outliers[sample_name][0]:
		pdb.set_trace()
	if line_outlier[1] != outliers[sample_name][1]:
		pdb.set_trace()
	t.write(str(switched) + '\t')
	t.write('\t'.join(np.asarray(line_outlier)) + '\t' + '\t'.join(np.asarray(line_score)) + '\t')
	t.write('\t'.join(np.asarray(outliers[sample_name])) + '\t' + '\t'.join(np.asarray(scores[sample_name])) + '\n')
f.close()
t.close()






