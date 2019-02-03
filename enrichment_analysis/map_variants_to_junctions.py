import numpy as np 
import os
import sys
import pdb


# Don't want a list of clusters where a bunch of them are repeats!
def remove_duplicates(old_string, add_on):
	if old_string == 'NULL':
		return add_on
	else:
		arr = old_string.split(',')
		if add_on in arr:
			return old_string
		else:
			return old_string + ',' + add_on	


# Mark on chromosome that the posi splice site and distance window around it map to this cluster_id
def fill_in_chromosome(chrom, posi, cluster_id, leftward_distance, rightward_distance):
	new_start = posi - leftward_distance
	if new_start < 0:
		new_start = 0
	new_end = posi + rightward_distance
	for pos in range(new_start, new_end + 1):
		chrom[pos] = remove_duplicates(chrom[pos], cluster_id)
	return chrom


def get_tissues(file_name):
	arr = []
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	return arr


def make_chromosome(chromosome, outlier_file, chrom_string, distance):
	f = open(outlier_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		jxn_info = data[0].split(':')
		if jxn_info[0] != chrom_string:
			continue
		start = int(jxn_info[1])
		end = int(jxn_info[2])
		jxn_id = jxn_info[0] + ':' + jxn_info[1] + ':' + jxn_info[2]
		# Mark on chromosome that the start splice site and distance window around it map to this cluster_id
		chromosome = fill_in_chromosome(chromosome, start, jxn_id, distance-1, distance)
		# Mark on chromosome that the end splice site and distance window around it map to this cluster_id
		chromosome = fill_in_chromosome(chromosome, end, jxn_id, distance, distance-1)
	f.close()
	return chromosome

def stream_variant_bed_file(chromosome, variant_bed_file, t,chrom_num):
	f = open(variant_bed_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[0] != chrom_num:
			continue
		pos = int(data[1])
		jxn_stringer = chromosome[pos]
		if jxn_stringer == 'NULL':
			continue
		jxns = np.unique(jxn_stringer.split(','))
		for jxn in jxns:
			t.write(line + '\t' + jxn + '\n')
	f.close()
	return t



variant_bed_file = sys.argv[1]
variant_junction_bed_file = sys.argv[2]
heuristic_outlier_dir = sys.argv[3]
heuristic_outlier_suffix = sys.argv[4]
tissue_names_file = sys.argv[5]
distance = int(sys.argv[6])

# Open output handle
t = open(variant_junction_bed_file,'w')

tissues = get_tissues(tissue_names_file)


for chrom_num in range(1,23):
	print(chrom_num)
	chromosome = ['NULL']*259250621
	for tissue in tissues:
		outlier_file = heuristic_outlier_dir + tissue + heuristic_outlier_suffix
		chromosome = make_chromosome(chromosome, outlier_file, 'chr' + str(chrom_num), distance)
	t = stream_variant_bed_file(chromosome, variant_bed_file, t, 'chr' + str(chrom_num))