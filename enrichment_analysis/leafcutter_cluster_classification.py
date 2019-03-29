import numpy as np 
import os
import sys
import pdb
from itertools import combinations



def there_is_exon_skipping_ordered(jxn_1, jxn_2, jxn_3):
	temper_bool = False
	jxn_1_start = int(jxn_1.split(':')[0])
	jxn_1_end = int(jxn_1.split(':')[1])
	jxn_2_start = int(jxn_2.split(':')[0])
	jxn_2_end = int(jxn_2.split(':')[1])
	jxn_3_start = int(jxn_3.split(':')[0])
	jxn_3_end = int(jxn_3.split(':')[1])
	if jxn_1_start == jxn_2_start and jxn_1_end == jxn_3_end and jxn_2_end < jxn_3_start:
		temper_bool = True
	if jxn_1_start == jxn_3_start and jxn_1_end == jxn_2_end and jxn_3_end < jxn_2_start:
		temper_bool = True
	return temper_bool

def there_is_alternative_5_prime(jxn_1, jxn_2):
	temper_bool = False
	jxn_1_start = int(jxn_1.split(':')[0])
	jxn_1_end = int(jxn_1.split(':')[1])
	jxn_2_start = int(jxn_2.split(':')[0])
	jxn_2_end = int(jxn_2.split(':')[1])
	if jxn_1_end == jxn_2_end and jxn_1_start != jxn_2_start:
		temper_bool = True
	return temper_bool

def there_is_alternative_3_prime(jxn_1, jxn_2):
	temper_bool = False
	jxn_1_start = int(jxn_1.split(':')[0])
	jxn_1_end = int(jxn_1.split(':')[1])
	jxn_2_start = int(jxn_2.split(':')[0])
	jxn_2_end = int(jxn_2.split(':')[1])
	if jxn_1_start == jxn_2_start and jxn_1_end != jxn_2_end:
		temper_bool = True
	return temper_bool

# Test whether this triplet of jxns exhibits exon skipping
def there_is_exon_skipping(jxn_1, jxn_2, jxn_3):
	temp_bool = False
	if there_is_exon_skipping_ordered(jxn_1, jxn_2, jxn_3):
		temp_bool = True
	if there_is_exon_skipping_ordered(jxn_2, jxn_1, jxn_3):
		temp_bool = True
	if there_is_exon_skipping_ordered(jxn_3, jxn_1, jxn_2):
		temp_bool = True
	return temp_bool

def check_if_jxn_pair_is_in_exon_skipping_triplet(jxn_1_index, jxn_2_index, exon_skipping_triplets):
	temp_bool = False
	for triplet_dicti in exon_skipping_triplets:
		if jxn_1_index in triplet_dicti and jxn_2_index in triplet_dicti:
			temp_bool = True
	return temp_bool


# Classify a leafcutter cluster into a specific splicing class
def classify_cluster(exon_exon_junctions, strand):
	# Initialize Binary variable explaining cluster is a certain classification
	exon_skipping = False
	alternative_5 = False
	alternative_3 = False
	# Need to loop through all N choose 3 triplets of junctions to identify if there is any exon skipping
	num_jxns = len(exon_exon_junctions)
	comb_3 = combinations(range(num_jxns), 3)
	exon_skipping_triplets = []
	for i in list(comb_3):
		jxn_1 = exon_exon_junctions[i[0]]
		jxn_2 = exon_exon_junctions[i[1]]
		jxn_3 = exon_exon_junctions[i[2]]
		if there_is_exon_skipping(jxn_1, jxn_2, jxn_3):
			exon_skipping = True
			exon_skipping_triplets.append({i[0]:1, i[1]:1, i[2]:1})

	# Need to loop through all N choose 2 pairs of junctions to identify if there is any alternate 5 and 3 sites
	comb_2 = combinations(range(num_jxns), 2)
	for i in list(comb_2):
		jxn_1 = exon_exon_junctions[i[0]]
		jxn_2 = exon_exon_junctions[i[1]]
		if check_if_jxn_pair_is_in_exon_skipping_triplet(i[0], i[1], exon_skipping_triplets):
			continue
		if there_is_alternative_5_prime(jxn_1, jxn_2):
			if strand == '+':
				alternative_5 = True
			elif strand == '-':
				alternative_3 = True
			else:
				print('assumption error!')
		if there_is_alternative_3_prime(jxn_1, jxn_2):
			if strand == '-':
				alternative_5 = True
			elif strand == '+':
				alternative_3 = True
			else:
				print('assumption error!')
	return exon_skipping, alternative_5, alternative_3