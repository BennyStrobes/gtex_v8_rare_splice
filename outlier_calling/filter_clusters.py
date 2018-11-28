import numpy as np 
import os
import sys
import pdb
import gzip


#produce list of autosomal chromosomes
def get_valid_chromosomes():
    valid_chromosomes = {}
    for num in range(1, 23):
        stringer = 'chr' + str(num)
        valid_chromosomes[stringer] = 1
    return valid_chromosomes


# Extract integer list of number of reads mapping to this jxn in each sample
def extract_jxn_reads(data):
	num_reads = []
	for ele in data:
		read_count = int(ele.split('/')[0])
		num_reads.append(read_count)
	return np.asarray(num_reads)

# Return true if one sample has at least min reads
def check_if_jxn_has_sample_with_more_than_min_reads(jxn_reads, min_reads):
	if max(jxn_reads) >= min_reads:
		return True 
	else:
		return False

# Create mapping from cluster id to a count of the number of junctions in the tissue that pass the filters
def create_mapping_from_cluster_id_to_num_jxns(raw_leafcutter_cluster_file, min_reads):
	# Create dictionary list of valid chromsome strings
	valid_chromosomes = get_valid_chromosomes()

	#mapping from cluster id to a count of the number of junctions in the tissue that pass the filters
	cluster_id_to_num_jxn = {}
	#mapping from cluster id to a vector (of length number of samples) that shows number of reads summed across all jxns for that cluster
	cluster_id_to_num_reads = {}

	f = gzip.open(raw_leafcutter_cluster_file)
	head_count = 0 # Used to skip header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Standard non-header line
		jxn_name = data[0]

		# Extract integer list of number of reads mapping to this jxn in each sample
		jxn_reads = extract_jxn_reads(data[1:])

		# Parse jxn name
		jxn_info = jxn_name.split(':')
		chromer = jxn_info[0]
		cluster_id = jxn_info[3]

		# Check if jxn passes our filters
		if chromer in valid_chromosomes and check_if_jxn_has_sample_with_more_than_min_reads(jxn_reads, min_reads):
			# Keep track of number of jxns that pass filters that this cluster has 
			if cluster_id not in cluster_id_to_num_jxn:
				cluster_id_to_num_jxn[cluster_id] = 1
				cluster_id_to_num_reads[cluster_id] = jxn_reads
			else:
				cluster_id_to_num_jxn[cluster_id] = cluster_id_to_num_jxn[cluster_id] + 1
				cluster_id_to_num_reads[cluster_id] = cluster_id_to_num_reads[cluster_id] + jxn_reads
	f.close()
	return cluster_id_to_num_jxn, cluster_id_to_num_reads


# Check if cluster has at least $min_samples_per_cluster with $min_reads_per_sample_in_cluster summed across all valid junctions
def check_if_cluster_has_enough_samples(num_reads, min_reads_per_sample_in_cluster, min_samples_per_cluster):
	if len(np.where(num_reads >= min_reads_per_sample_in_cluster)[0]) >= min_samples_per_cluster:
		return True
	else:
		return False

def filter_jxn_file(raw_leafcutter_cluster_file, filtered_leafcutter_cluster_file, min_reads, min_reads_per_sample_in_cluster, min_samples_per_cluster, cluster_id_to_num_jxn, cluster_id_to_num_reads):
	# Create dictionary list of valid chromsome strings
	valid_chromosomes = get_valid_chromosomes()
	# Open output file
	t = open(filtered_leafcutter_cluster_file, 'w')
	# Open input file
	f = gzip.open(raw_leafcutter_cluster_file)
	head_count = 0 # Used to skip header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			# Print header (though with individual ids instead of gtex sample ids)
			t.write(data[0])
			for ele in data[1:]:
				individual_id = ele.split('-')[0] + '-' + ele.split('-')[1]
				t.write('\t' + individual_id)
			t.write('\n')
			continue
		# Standard non-header line
		jxn_name = data[0]

		# Extract integer list of number of reads mapping to this jxn in each sample
		jxn_reads = extract_jxn_reads(data[1:])

		# Parse jxn name
		jxn_info = jxn_name.split(':')
		chromer = jxn_info[0]
		cluster_id = jxn_info[3]

		# Check if jxn passes our filters
		if chromer in valid_chromosomes and check_if_jxn_has_sample_with_more_than_min_reads(jxn_reads, min_reads):
			# Check if cluster passes our filters (cluster has at least 2 jxns that pass filters) AND at least $min_samples_per_cluster with $min_reads_per_sample_in_cluster summed across all valid junctions
			if cluster_id_to_num_jxn[cluster_id] >= 2 and check_if_cluster_has_enough_samples(cluster_id_to_num_reads[cluster_id], min_reads_per_sample_in_cluster, min_samples_per_cluster):
				t.write(jxn_name + '\t' + '\t'.join(jxn_reads.astype(str)) + '\n')
	f.close()
	t.close()


raw_leafcutter_cluster_file = sys.argv[1]
filtered_leafcutter_cluster_file = sys.argv[2]
min_reads = int(sys.argv[3])
min_reads_per_sample_in_cluster = int(sys.argv[4])
min_samples_per_cluster = int(sys.argv[5])

# Create mapping from cluster id to a count of the number of junctions in the tissue that pass the filters (autosomal chromosome and at least 1 sample has at least $min_reads)
# Create mapping from cluster id to a vector (of length number of samples) that shows number of reads summed across all valid jxns for that cluster
cluster_id_to_num_jxn, cluster_id_to_num_reads = create_mapping_from_cluster_id_to_num_jxns(raw_leafcutter_cluster_file, min_reads)

# Filter $raw_leafcutter_cluster_file and print results to $filtered_leafcutter_cluster_file
filter_jxn_file(raw_leafcutter_cluster_file, filtered_leafcutter_cluster_file, min_reads, min_reads_per_sample_in_cluster, min_samples_per_cluster, cluster_id_to_num_jxn, cluster_id_to_num_reads)
