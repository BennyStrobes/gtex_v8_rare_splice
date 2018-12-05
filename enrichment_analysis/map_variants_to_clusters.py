import numpy as np 
import os
import sys
import pdb

def quick_error_check_to_make_sure_there_are_at_least_2_jxns(jxns):
	if len(jxns) < 2:
		print('jxn assumption error')
		pdb.set_trace()
	return

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
def fill_in_chromosome(chrom, posi, cluster_id, distance):
	new_start = posi - distance
	if new_start < 0:
		new_start = 0
	new_end = posi + distance
	for pos in range(new_start, new_end + 1):
		chrom[pos] = remove_duplicates(chrom[pos], cluster_id)
	return chrom

# Fill in the object for the current chromosome
# Keep track of which BP on the chromosome correspond to which cluster
def make_chromosome(cluster_info_file, chrom_num, distance):
	# Initialize
	chrom = ['NULL']*259250621
	chrom_string = 'chr' + str(chrom_num)
	# Stream cluster file
	f = open(cluster_info_file)
	head_count = 0  # For header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Standard line
		# Parse line
		cluster_id = data[0]
		jxns = data[1].split(',')
		chromer = jxns[0].split(':')[0]
		# Skip clusters not on this chromosome
		if chromer != chrom_string:
			continue
		quick_error_check_to_make_sure_there_are_at_least_2_jxns(jxns)
		# Loop through jxns for this cluster
		for jxn in jxns:
			# Parse jxn
			jxn_info = jxn.split(':')
			start = int(jxn_info[1])
			end = int(jxn_info[2])
			# Mark on chromosome that the start splice site and distance window around it map to this cluster_id
			chrom = fill_in_chromosome(chrom, start, cluster_id, distance)
			# Mark on chromosome that the end splice site and distance window around it map to this cluster_id
			chrom = fill_in_chromosome(chrom, end, cluster_id, distance)
	f.close()
	return chrom


# Stream variant bed file. For each variant on the current chromosome, try to map to a cluster
def stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, t_filter, chrom_num):
	chrom_string = 'chr' + str(chrom_num)
	# Stream variant file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		line_chrom_num = data[1]
		# Throw out variants not on current chromosome
		if line_chrom_num != chrom_string:
			continue
		# Make sure MAF is not >= .01
		maf = float(data[4])
		if maf >= .01:
			print('maf assumption error!')
			continue
		var_pos = int(data[2])  # Position of variant
		# Get list of cluster_ids that overlap the variant position
		overlapping_clusters = cluster_chromosome[var_pos].split(',')
		# Print one line for every variant-cluster pair
		# If no clusters map, put one line and a NULL where cluster_id goes
		for cluster_id in overlapping_clusters:
			t.write(line + '\t' + cluster_id + '\n')
			if cluster_id != 'NULL':
				t_filter.write(line + '\t' + cluster_id + '\n')
	f.close()
	return t, t_filter


##################################
# Command Line ARGS
##################################
variant_bed_file = sys.argv[1]  # Input variant file
variant_cluster_bed_file = sys.argv[2]  # Output variant file (will include info on mapping of variant to clusters)
variant_cluster_only_bed_file = sys.argv[3]  # output variant file (same as above but filters out variants not mapped to a cluster)
cluster_info_file = sys.argv[4]  # File containing info/location of each cluster
distance = int(sys.argv[5])  # Distance window around splice site that we want to consider variants


# Open output file handle
t = open(variant_cluster_bed_file, 'w')
t_filter = open(variant_cluster_only_bed_file, 'w')

# Loop through autosomal chromsosomes
for chrom_num in range(1,23):
	# Initialize object to keep track of which BP on the chromosome correspond to which cluster
	cluster_chromosome = {}
	# Fill in the object for the current chromosome
	cluster_chromosome = make_chromosome(cluster_info_file, chrom_num, distance)
	# Stream variant bed file. For each variant on the current chromosome, try to map to a cluster
	t, t_filter = stream_variant_bed_file(cluster_chromosome, variant_bed_file, t, t_filter, chrom_num)

#  Close file-handles
t.close()
t_filter.close()