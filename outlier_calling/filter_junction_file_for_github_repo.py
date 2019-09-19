import numpy as np 
import os
import sys
import pdb



def randomly_select_n_clusters(outlier_pvalue_file, num_clusters):
	clusters = {}
	f = open(outlier_pvalue_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			total_indi = len(data) -1
			continue
		cluster_id = data[0]
		pvalues = np.asarray(data[1:]).astype(float)
		if len(np.where(pvalues < .01)[0]) < 6:
			clusters[cluster_id] = []
	np.random.seed(1)
	n_clusters = np.random.choice(clusters.keys(), size=num_clusters, replace=False)
	dicti = {}
	for cluster in n_clusters:
		dicti[cluster] = []
	return dicti, total_indi


jxn_file = sys.argv[1]
outlier_pvalue_file = sys.argv[2]
filtered_jxn_file = sys.argv[3]


num_clusters = 20
num_indi = 300
clusters, total_indi = randomly_select_n_clusters(outlier_pvalue_file, num_clusters)

cluster_samples = {}
for cluster_id in clusters.keys():
	cluster_samples[cluster_id] = np.random.choice(range(total_indi), replace=False, size=num_indi)

f = open(jxn_file)
t = open(filtered_jxn_file,'w')

header = []
for indi_num in range(num_indi):
	header.append('sample_' + str(indi_num))

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(data[0] + '\t' + '\t'.join(header) + '\n')
		continue
	jxn_id = data[0]
	cluster_id = jxn_id.split(':')[3]
	if cluster_id not in clusters:
		continue
	counts = data[1:]
	subset_counts = np.asarray(counts)[cluster_samples[cluster_id]]
	new_line = data[0] + '\t' + '\t'.join(subset_counts)
	clusters[cluster_id].append(new_line)
f.close()
for cluster in clusters.keys():
	for liner in clusters[cluster]:
		t.write(liner + '\n')
t.close()