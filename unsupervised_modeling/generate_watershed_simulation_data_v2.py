import numpy as np
import os
import sys
import pdb
import scipy.stats



def generate_new_pvalues(pvalue):
	if pvalue > .1:
		new_pvalues = np.random.uniform(low=.1,high=1,size=3)
	else:
		randy =np.random.uniform(low=0,high=1,size=1)
		if randy > .6:
			new_pvalues = np.random.uniform(low=0,high=.1,size=3)
		elif randy > .3:
			new_pvalues = np.random.uniform(low=0,high=.1,size=3)
			new_pvalues[0] = np.random.uniform(low=0,high=1,size=1)
		elif randy > .2:
			new_pvalues = np.random.uniform(low=0,high=.1,size=3)
			new_pvalues[1] = np.random.uniform(low=0,high=1,size=1)
		else:
			new_pvalues = np.random.uniform(low=0,high=.3,size=3)
	return new_pvalues

def add_specific_feats(feat_vec, new_pvalues, feat_1, feat_2, feat_3):
	new_vec = np.asarray(feat_vec).astype(float)
	if new_pvalues[0] < .1 and np.random.uniform(low=0,high=1,size=1) < .01:
		new_vec = feat_1 + np.random.normal(scale=.6,size=18) 
	if new_pvalues[1] < .1 and np.random.uniform(low=0,high=1,size=1) < .01:
		new_vec = feat_2 + np.random.normal(scale=.6,size=18)
	if new_pvalues[2] < .1 and np.random.uniform(low=0,high=1,size=1) < .01:
		new_vec = feat_3 + np.random.normal(scale=.6,size=18)
	return np.asarray(new_vec).astype(str)

def reduce_watershed_to_river(watershed_sim_file, pheno_num, river_sim_file_stem):
	output_file = river_sim_file_stem + 'pheno_' + str(pheno_num) + '.txt'
	f = open(watershed_sim_file)
	t = open(output_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write('\t'.join(data[0:20]) + '\t')
		t.write(data[(19+pheno_num)] + '\t' + data[23] + '\n')
	f.close()
	t.close()


river_sim_file = sys.argv[1]
watershed_sim_file = sys.argv[2]
river_sim_file_stem = sys.argv[3]
'''
np.random.seed(8)
num_feat = 18
feat_1 = np.random.normal(scale=.1,size=18)
feat_2 = np.random.normal(scale=.1,size=18)
feat_3 = np.random.normal(scale=.1,size=18)

f = open(river_sim_file)
t = open(watershed_sim_file, 'w')
dicti = {}
dist = []
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write('\t'.join(data[0:20]) + '\tpheno_1_pvalue\tpheno_2_pvalue\tpheno_3_pvalue\tN2pair\n')
		continue
	zscore = float(data[-2])
	pvalue = scipy.stats.norm.sf(abs(zscore))*2
	new_pvalues = generate_new_pvalues(pvalue)
	new_feat = add_specific_feats(data[2:20], new_pvalues, feat_1, feat_2, feat_3)
	#new_feat = np.asarray(data[2:20]).astype(str)
	if zscore < 0:
		new_pvalues[1] = new_pvalues[1]*-1
	if data[-1] != 'NA':
		if data[-1] not in dicti:
			dicti[data[-1]] = (new_pvalues, new_feat)
		else:
			new_feat = dicti[data[-1]][1]
			if float(dicti[data[-1]][0][0]) < .1 and abs(float(dicti[data[-1]][0][1])) < .1 and float(dicti[data[-1]][0][2]) < .1:
				randy =np.random.uniform(low=0,high=1,size=1)
				new_pvalues = np.random.uniform(low=0,high=.1,size=3)
			elif abs(float(dicti[data[-1]][0][1])) < .1 and float(dicti[data[-1]][0][2]) < .1:
				randy =np.random.uniform(low=0,high=1,size=1)
				new_pvalues = np.random.uniform(low=0,high=.1,size=3)
				new_pvalues[0] = np.random.uniform(low=0,high=1,size=1)
			elif float(dicti[data[-1]][0][0]) < .1 and float(dicti[data[-1]][0][2]) < .1:
				randy =np.random.uniform(low=0,high=1,size=1)
				new_pvalues = np.random.uniform(low=0,high=.1,size=3)
				new_pvalues[1] = np.random.uniform(low=0,high=1,size=1)
			elif float(dicti[data[-1]][0][0]) < .1 and abs(float(dicti[data[-1]][0][1])) < .1:
				new_pvalues = np.random.uniform(low=0,high=.1,size=3)
				new_pvalues[2] = np.random.uniform(low=0,high=1,size=1)
	t.write(data[0] + '\t' + data[1] + '\t' + '\t'.join(new_feat) + '\t' + '\t'.join(new_pvalues.astype(str)) + '\t' + data[21] + '\n')


t.close()
'''

reduce_watershed_to_river(watershed_sim_file, 1, river_sim_file_stem)
reduce_watershed_to_river(watershed_sim_file, 2, river_sim_file_stem)
reduce_watershed_to_river(watershed_sim_file, 3, river_sim_file_stem)

