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

# Create mapping from cluster id to a count of the number of junctions in that cluster
# Create mapping from cluster id to a vector (of length number of samples) that shows number of reads summed across all valid jxns for that cluster
def create_mapping_from_cluster_id_to_num_jxns(raw_leafcutter_cluster_file, min_reads):
	# Create dictionary list of valid chromsome strings
	valid_chromosomes = get_valid_chromosomes()

	#mapping from cluster id to a count of the number of junctions in the tissue that pass the filters
	cluster_id_to_num_jxn = {}
	#mapping from cluster id to a vector (of length number of samples) that shows number of reads summed across all jxns for that cluster
	cluster_id_to_num_reads = {}

	f = open(raw_leafcutter_cluster_file)
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

# Filter $raw_leafcutter_cluster_file and print results to $filtered_leafcutter_cluster_file
# Remove clusters that:
## 1. Only have 1 junction
## 2. Have fewer than $min_samples_per_cluster with >= $min_reads_per_sample_in_cluster
def filter_jxn_file(raw_leafcutter_cluster_file, filtered_leafcutter_cluster_file, min_reads, min_reads_per_sample_in_cluster, min_samples_per_cluster, cluster_id_to_num_jxn, cluster_id_to_num_reads):
	# Create dictionary list of valid chromsome strings
	valid_chromosomes = get_valid_chromosomes()
	# Open output file
	t = open(filtered_leafcutter_cluster_file, 'w')
	# Open input file
	f = open(raw_leafcutter_cluster_file)
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

# The cluster ss1 is previously mapped to and the cluster ss2 is previously mapped to disagree. --> Merge clusters to a new one and delete old
def merge_clusters(ss1, ss2, ss_to_clusters, clusters, cluster_number, header_line):
    cluster_name_old_1 = ss_to_clusters[ss1]  # Cluster name corresponding to old ss1
    cluster_name_old_2 = ss_to_clusters[ss2]  # Cluster name corresponding to old ss2
    cluster_name = 'cluster' + str(cluster_number)  # New cluster_id

    ss_to_clusters[ss1] = cluster_name  # Map ss1 to new_cluster_id
    ss_to_clusters[ss2] = cluster_name  # Map ss2 to new_cluster_id
    clusters[cluster_name] = []
    clusters[cluster_name].append(header_line)  # Map cluster_id to this jxn

    for old_header_line in clusters[cluster_name_old_1]:  # Remap all jxns in old cluster name corresponding to ss1 to new cluster id
        line_header_info = old_header_line.split(':')
        line_ss1 = line_header_info[0] + '_' + line_header_info[1]
        line_ss2 = line_header_info[0] + '_' + line_header_info[2]
        ss_to_clusters[line_ss1] = cluster_name
        ss_to_clusters[line_ss2] = cluster_name
        clusters[cluster_name].append(old_header_line)
    for old_header_line in clusters[cluster_name_old_2]:  # Remap all jxns in old cluster name corresponding to ss2 to new cluster id
        line_header_info = old_header_line.split(':')
        line_ss1 = line_header_info[0] + '_' + line_header_info[1]
        line_ss2 = line_header_info[0] + '_' + line_header_info[2]
        ss_to_clusters[line_ss1] = cluster_name
        ss_to_clusters[line_ss2] = cluster_name
        clusters[cluster_name].append(old_header_line)

    clusters.pop(cluster_name_old_1)  # Remove old cluster names
    clusters.pop(cluster_name_old_2)  # Remove old cluster names
    return ss_to_clusters, clusters

def check_to_make_sure_individuals_match_other_rv_individauls(sample_names, individual_list, tissue_name):
	known_individuals = {}
	f = open(individual_list)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		tiss_name = data[0]
		indi_id = data[1]
		if tiss_name != tissue_name:
			continue
		known_individuals[indi_id] = 0
	f.close()
	for sample in sample_names:
		indi = sample.split('-')[0] + '-' + sample.split('-')[1]
		if indi not in known_individuals:
			print('erroror!!')
			pdb.set_trace()
		known_individuals[indi] = 1
	for indi in known_individuals:
		if known_individuals[indi] != 1:
			print('erororor')
			pdb.set_trace()
	return

# Remove Junctions that are non-autosomal and don't have at least 1 sample with at least $min_reads
def filter_junctions_and_re_assign_clusters(raw_leafcutter_cluster_file, temp_raw_cluster_file, min_reads, individual_list, tissue_name):
	# Mapping from cluster ID to jxns in that cluster
	clusters = {}
	# Mapping from SS to cluster ID
	ss_to_clusters = {}
	# Integer keeping track of current cluster
	cluster_number = 0

	# Create dictionary list of valid chromsome strings
	valid_chromosomes = get_valid_chromosomes()
	# Open input file
	f = gzip.open(raw_leafcutter_cluster_file)
	head_count = 0 # Used to skip header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			check_to_make_sure_individuals_match_other_rv_individauls(data[1:], individual_list, tissue_name)
			continue
		# Standard non-header line
		jxn_name = data[0]

		# Extract integer list of number of reads mapping to this jxn in each sample
		jxn_reads = extract_jxn_reads(data[1:])

		# Parse jxn name
		jxn_info = jxn_name.split(':')
		chromer = jxn_info[0]
		cluster_id = jxn_info[3]
		ss1 = jxn_info[0] + '_' + jxn_info[1]
		ss2 = jxn_info[0] + '_' +jxn_info[2]
		# Quick error check
		if int(jxn_info[1]) > int(jxn_info[2]):
			print('SS assumption error')
			pdb.set_trace()

		# Ignore Jxns not on valid chromosomes
		if chromer not in valid_chromosomes:
			continue
		# Ignore jxns with no samples with >= min_reads
		if check_if_jxn_has_sample_with_more_than_min_reads(jxn_reads, min_reads) == False:
			continue

		# Asssign jxn to a cluster
		# Four different cases to deal with...
		if ss1 not in ss_to_clusters and ss2 not in ss_to_clusters:  # Neither splice site has been seen before --> create new cluster_id
			cluster_name = 'cluster' + str(cluster_number)
			ss_to_clusters[ss1] = cluster_name
			ss_to_clusters[ss2] = cluster_name
			clusters[cluster_name] = []
			clusters[cluster_name].append(jxn_name)
			cluster_number = cluster_number + 1
		elif ss1 not in ss_to_clusters and ss2 in ss_to_clusters:  # ss2 has been seen before, but ss1 has not been seen before. Map this jxn to the cluster ss2 is mapped to
			cluster_name = ss_to_clusters[ss2]
			ss_to_clusters[ss1] = cluster_name
			ss_to_clusters[ss2] = cluster_name
			clusters[cluster_name].append(jxn_name)
		elif ss1 in ss_to_clusters and ss2 not in ss_to_clusters:  # ss1 has been seen before, but ss2 has not been seen before. Map this jxn to the cluster ss1 is mapped to.
			cluster_name = ss_to_clusters[ss1]
			ss_to_clusters[ss1] = cluster_name
			ss_to_clusters[ss2] = cluster_name
			clusters[cluster_name].append(jxn_name)
		elif ss1 in ss_to_clusters and ss2 in ss_to_clusters:  # ss1 and ss2 have been seen before (most interesting case)
			cluster_name1 = ss_to_clusters[ss1]
			cluster_name2 = ss_to_clusters[ss2]
			if cluster_name1 != cluster_name2:  # The cluster ss1 is previously mapped to and the cluster ss2 is previously mapped to disagree. --> Merge clusters to a new one and delete old
				ss_to_clusters, clusters = merge_clusters(ss1, ss2, ss_to_clusters, clusters, cluster_number, jxn_name)
				cluster_number = cluster_number + 1
			else:  # ss1 and ss2 previously mapped to the same cluster.
				cluster_name = cluster_name1
				ss_to_clusters[ss1] = cluster_name
				ss_to_clusters[ss2] = cluster_name
				clusters[cluster_name].append(jxn_name)
	f.close()

	# Now that we have mapped valid junctions to new cluster assignments, print

	# Open input file
	f = gzip.open(raw_leafcutter_cluster_file)
	# Open output file
	t = open(temp_raw_cluster_file, 'w')
	head_count = 0 # Used to skip header
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		# Standard non-header line
		jxn_name = data[0]

		# Extract integer list of number of reads mapping to this jxn in each sample
		jxn_reads = extract_jxn_reads(data[1:])

		# Parse jxn name
		jxn_info = jxn_name.split(':')
		chromer = jxn_info[0]
		cluster_id = jxn_info[3]

		# Ignore Jxns not on valid chromosomes
		if chromer not in valid_chromosomes:
			continue
		# Ignore jxns with no samples with >= min_reads
		if check_if_jxn_has_sample_with_more_than_min_reads(jxn_reads, min_reads) == False:
			continue

		# Splice site names
		ss1 = jxn_info[0] + '_' + jxn_info[1]
		ss2 = jxn_info[0] + '_' +jxn_info[2]

		if ss1 not in ss_to_clusters or ss2 not in ss_to_clusters:
			pdb.set_trace()
		# Map from SS to new cluster assignment
		cluster_name1 = ss_to_clusters[ss1]
		cluster_name2 = ss_to_clusters[ss2]
		if cluster_name1 != cluster_name2:
			print('EROROOROR')
			pdb.set_trace()
		new_header_line = jxn_info[0] + ':' + jxn_info[1] + ':' + jxn_info[2] + ':' + cluster_name1  # Jxnid string for output
		t.write(new_header_line + '\t' + '\t'.join(data[1:]) + '\n')
	f.close()
	t.close()



raw_leafcutter_cluster_file = sys.argv[1]
filtered_leafcutter_cluster_file = sys.argv[2]
min_reads = int(sys.argv[3])
min_reads_per_sample_in_cluster = int(sys.argv[4])
min_samples_per_cluster = int(sys.argv[5])
individual_list = sys.argv[6]
tissue_name = sys.argv[7]

# Intermediate junction file that will be deleted at the end of this script
temp_raw_cluster_file = filtered_leafcutter_cluster_file.split('.')[0] + '_temp.txt' 

# Remove Junctions that are non-autosomal and don't have at least 1 sample with at least $min_reads
# Also check to make sure list of individuals being used matches other RV analysis individual lists
filter_junctions_and_re_assign_clusters(raw_leafcutter_cluster_file, temp_raw_cluster_file, min_reads, individual_list, tissue_name)

# Create mapping from cluster id to a count of the number of junctions in that cluster
# Create mapping from cluster id to a vector (of length number of samples) that shows number of reads summed across all valid jxns for that cluster
cluster_id_to_num_jxn, cluster_id_to_num_reads = create_mapping_from_cluster_id_to_num_jxns(temp_raw_cluster_file, min_reads)

# Filter $raw_leafcutter_cluster_file and print results to $filtered_leafcutter_cluster_file
# Remove clusters that:
## 1. Only have 1 junction
## 2. Have fewer than $min_samples_per_cluster with >= $min_reads_per_sample_in_cluster
filter_jxn_file(temp_raw_cluster_file, filtered_leafcutter_cluster_file, min_reads, min_reads_per_sample_in_cluster, min_samples_per_cluster, cluster_id_to_num_jxn, cluster_id_to_num_reads)


# Delete intermediate junction file
os.system('rm ' + temp_raw_cluster_file)
