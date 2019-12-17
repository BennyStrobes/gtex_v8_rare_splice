import numpy as np 
import os
import gzip
import sys
import pdb


def get_cluster_to_strand_mapping(cluster_info_file, exon_file):
	gene_to_strand = {}
	cluster_to_strand = {}
	f = open(exon_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ensamble_id = data[4]
		strand = data[3]
		gene_to_strand[ensamble_id] = strand
	f.close()
	head_count = 0
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		ensamble_ids = data[2].split(',')
		for ensamble_id in ensamble_ids:
			cluster_to_strand[cluster_id] = gene_to_strand[ensamble_id]
	f.close()
	return cluster_to_strand

def correct_alleles_for_strand(strand, major_allele, variant_allele):
	if strand == '-':
		converter = {'A':'T', 'T':'A','G':'C','C':'G'}
		return converter[major_allele], converter[variant_allele]
	else:
		return major_allele, variant_allele
def get_to_or_from_concensus(concensus_allele, variant_allele):
	old = variant_allele.split('->')[0]
	new = variant_allele.split('->')[1]
	if old == concensus_allele:
		return 'from_concensus'
	elif new == concensus_allele:
		return 'to_concensus'
	else:
		return 'neither'

def run_step_1_filters(inlier_variant_position_file, cluster_to_strand_mapping, step_1_output_file):
	valid_positions = {'A-2':'A','A-1':'G', 'A+1':'G', 'D-1':'G', 'D+1':'G', 'D+2':'T', 'D+3':'A', 'D+4':'A', 'D+5':'G', 'D+6':'T'}
	f = open(inlier_variant_position_file)
	t = open(step_1_output_file, 'w')
	head_count = 0
	t.write('variant_name\tindi_id\tcluster_id\tvariant_position\tto_consensus\n')
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		relative_position = float(data[0])
		if relative_position <= 0:
			continue
		if relative_position > 60:
			continue
		if data[1] == 'donor':
			variant_position = 'D-' + str(int(relative_position+1))
		else:
			variant_position = 'A+' + str(int(relative_position+1))
		variant_name = data[6] + '_' + data[7]
		cluster_id = data[11]
		indi_id = data[5]
		strand = cluster_to_strand_mapping[cluster_id]
		major_allele = data[3]
		variant_allele = data[4]
		major_allele, variant_allele = correct_alleles_for_strand(strand, major_allele, variant_allele)
		#variant_type = get_to_or_from_concensus(valid_positions[variant_position], major_allele + '->' + variant_allele)
		variant_type = major_allele + '->' + variant_allele
		t.write(variant_name + '\t' + indi_id + '\t' + cluster_id + '\t' + variant_position + '\t' + variant_type + '\n')
	f.close()
	t.close()

# Simple check to ensure everything that passed step 1's filters also is an outlier in  whole blood (should be!)
def run_step_2_filters(input_file, whole_blood_outlier_file, output_file, column_name):
	# First extract cluster_individual pairs that pass step 1 filter
	samples = {}
	clusters = {}
	indiz = {}
	head_count = 0
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		samples[data[2] + '_' + data[1]] = -1
		clusters[data[2]] = 0
	f.close()
	f = open(whole_blood_outlier_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			indis = np.asarray(data[1:])
			continue
		cluster_id = data[0]
		if cluster_id not in clusters:
			continue
		clusters[cluster_id] = 1
		pvalz = np.asarray(data[1:]).astype(float)
		for i, indi in enumerate(indis):
			pval = pvalz[i]
			if cluster_id + '_' + indi in samples:
				samples[cluster_id + '_' + indi] = pval
	f.close()
	t = open(output_file, 'w')
	f = open(input_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + column_name + '\n')
			continue
		pval = samples[data[2] + '_' + data[1]]
		if pval == -1:
			continue
		if pval > .35:
			t.write(line + '\t' + str(pval) + '\n')
	f.close()
	t.close()

def run_step_4_add_cluster_info(step_3_output_file, cluster_info_file, step_4_output_file):
	cluster_dicti = {}
	f = open(cluster_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		jxn_names = data[1]
		gene_names = data[2]
		cluster_dicti[cluster_id] = (jxn_names, gene_names)
	f.close()
	t = open(step_4_output_file, 'w')
	f = open(step_3_output_file)
	head_count= 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'jxn_names\tgene_names\n')
			continue
		cluster_id = data[2]
		info = cluster_dicti[cluster_id]
		if len(info[1].split(',')) == 1:
			t.write(line + '\t' + info[0] + '\t' + info[1] + '\n')
		else:
			print('skip cause multiple genes')
	f.close()
	t.close()

def run_step_5_add_watershed_info(step_4_output_file, watershed_posterior_file, step_5_output_file):
	samples = {}
	head_count = 0
	f = open(step_4_output_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_name = data[1] + ':' + data[8] + ':' + data[0]
		samples[sample_name] = 'hey'
	f.close()
	f = gzip.open(watershed_posterior_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count +1
			continue
		sample_name = data[0].split('_')[0] + '_' + data[0].split('_')[1]
		if sample_name in samples:
			samples[sample_name] = data
	f.close()
	f = open(step_4_output_file)
	t = open(step_5_output_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'GAM_splice\twatershed_splice\twatershed_ase\twatershed_expr\n')
			continue
		sample_name = data[1] + ':' + data[8] + ':' + data[0]
		watershed_info = samples[sample_name]
		if len(watershed_info[0].split(':')) != 3:
			continue
		variant_id = watershed_info[0].split(':')[2]
		# Print existing line (swapping in more informative variant id)
		if float(watershed_info[10]) < .2 and float(watershed_info[12]) < .2 and float(watershed_info[11]) < .2:
			t.write(variant_id + '\t' + '\t'.join(data[1:]))
			t.write('\t' + watershed_info[4] + '\t' + watershed_info[10] + '\t' + watershed_info[12] + '\t' + watershed_info[11] + '\n')
	t.close()
	f.close()

def run_step_6_check_to_see_if_variant_has_high_watershed_in_any_individuals(step_6_output_file, watershed_posterior_file, step_7_output_file):
	variants = {}
	head_count = 0
	variants = {}
	f = open(step_6_output_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[0] + ':' + data[8]
		variants[variant_id] = float("NaN")
	f.close()
	f = gzip.open(watershed_posterior_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count +1
			continue
		variant_id = data[0].split(':')[2] + ':' + data[0].split(':')[1]
		if variant_id in variants:
			watershed_splice = float(data[10])
			variants[variant_id] = np.nanmax((float(data[10]), variants[variant_id]))
	f.close()
	f = open(step_6_output_file)
	t = open(step_7_output_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'max_splicing_watershed_for_variant\n')
			continue
		variant_id = data[0] + ':' + data[8]
		if variants[variant_id] < .2:
			t.write(line + '\t' + str(variants[variant_id]) + '\n')
	f.close()
	t.close()

def organized_test_variants(step_7_output_file, organized_output_file, organized_with_junctions_output_file):
	f = open(step_7_output_file)
	t = open(organized_output_file, 'w')
	t2 = open(organized_with_junctions_output_file,'w')
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[8] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + '\t'.join(data[9:]) + '\n')
		t2.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + data[8] + '\t' + data[3] + '\t' + data[4] + '\t' + data[5] + '\t' + data[6] + '\t' + '\t'.join(data[9:]) + '\t' + data[7] + '\n')
	f.close()
	t.close()

# Create mapping from cluster id to junctions
def extract_mapping_from_cluster_id_to_junctions(cluster_info_file):
	# Initialize mapping
	mapping = {}
	# Used to skip header
	head_count = 0
	# Stream cluster_info_file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent info
		cluster_id = data[0]
		jxn_string = data[1]
		chromosome = jxn_string.split(',')[0].split(':')[0]
		jxn_array = jxn_string.split(',')
		# Add info to mapping
		if cluster_id in mapping:
			print('assumption error')
		mapping[cluster_id] = {}
		mapping[cluster_id]['chromosome'] = chromosome
		mapping[cluster_id]['start_splice_sites'] = []
		mapping[cluster_id]['end_splice_sites'] = []
		# Add splice sites to mapping
		for jxn in jxn_array:
			jxn_info = jxn.split(':')
			start_pos = int(jxn_info[1])
			end_pos = int(jxn_info[2])
			if end_pos <= start_pos:
				print('assumptioner eroror')
				pdb.set_trace()
			mapping[cluster_id]['start_splice_sites'].append(start_pos)
			mapping[cluster_id]['end_splice_sites'].append(end_pos)
	f.close()
	return mapping


# Create dictionary of annotated splice sites
# Where keys are "chr"$chrom_num"_"#splice_sites
def get_dictionary_of_annotated_splice_sites(exon_file):
	# initialize dictionary
	ss = {}
	# used to skip header
	head_count = 0
	# stream exon file
	f = open(exon_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		ss[data[0] + '_' + data[1]] = 1
		ss[data[0] + '_' + data[2]] = 1
	return ss

def get_inlier_indi(outlier_indi, all_individuals):
	inlier_indi = []
	count = 0
	for indi in all_individuals:
		if indi not in outlier_indi:
			inlier_indi.append(indi)
		else:
			count = count + 1
	if count != len(outlier_indi):
		print('assumption error')
		pdb.set_trace()
	return np.asarray(inlier_indi)

def make_clusters_to_plot_file(organized_output_file, cluster_to_ss_mapping, cluster_to_strand_mapping, annotated_splice_sites, whole_blood_outlier_file, clusters_to_plot_file):
	clusters = {}
	f = open(organized_output_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[2]
		variant_id = data[0]
		splice_site = data[0].split('_')[0] + ':' + data[0].split('_')[1]
		if cluster_id + '_' + variant_id not in clusters:
			if splice_site.split(':')[0] + '_' + splice_site.split(':')[1] in annotated_splice_sites:
				anno = 'annotated'
			else:
				anno = 'novel'
			tupler = ([data[1]], data[4], data[5], cluster_to_strand_mapping[cluster_id], anno, data[0])
			clusters[cluster_id + '_' + variant_id] = tupler
		else:
			if splice_site.split(':')[0] + '_' + splice_site.split(':')[1] in annotated_splice_sites:
				anno = 'annotated'
			else:
				anno = 'novel'
			old_tupler = clusters[cluster_id + '_' + variant_id]
			old_indi = old_tupler[0]
			old_indi.append(data[1])
			if data[4] != old_tupler[1]:
				print('assumption error')
				pdb.set_trace()
			if data[5] != old_tupler[2]:
				print('assumption error')
				pdb.set_trace()
			if cluster_to_strand_mapping[cluster_id] != old_tupler[3]:
				print('assumption eororo')
				pdb.set_trace()
			if anno != old_tupler[4]:
				print('assumption eror')
				pdb.set_trace()
			if data[0] != old_tupler[5]:
				print('assumption error')
				pdb.set_trace()
			tupler = (old_indi, data[4], data[5], cluster_to_strand_mapping[cluster_id], anno, data[0])
			clusters[cluster_id + '_' + variant_id] = tupler
	f.close()
	# Extract array of whole blood individuals
	all_individuals = []
	head_count = 0
	f = open(whole_blood_outlier_file)
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			head_count = head_count + 1
			all_individuals = np.asarray(data[1:])
			continue
	f.close()
	# print clusters_to_plot_file
	t = open(clusters_to_plot_file, 'w')
	t.write('cluster_id\tchrom_num\tvariant_position\toutlier_individual\tinlier_individuals\tdistance_to_ss\tstrand\tannotated_ss\tvariant_allele\n')
	# loop throug clusters
	# AFTER MANUAL FILTERING WE FOUND THIS SET OF CLUSTERS TO USE
	valid_clusters = {}
	valid_clusters['cluster14314_chr7_100367111_C_T_1'] = 1
	valid_clusters['cluster44749_chr14_92121331_A_C_1'] = 1
	valid_clusters['cluster56169_chr11_67277707_G_A_1'] = 1
	valid_clusters['cluster17388_chr4_2469184_T_C_1'] = 1
	valid_clusters['cluster34796_chr17_78113582_T_C_1'] = 1
	for cluster_id in clusters.keys():
		cluster_info = clusters[cluster_id]
		outlier_indi = cluster_info[0]
		distance_to_ss = cluster_info[1]
		variant_allele = cluster_info[2]
		strand = cluster_info[3]
		anno = cluster_info[4]
		variant_id = cluster_info[5]
		chrom_num = variant_id.split('_')[0]
		variant_position = variant_id.split('_')[1]
		inlier_indi = get_inlier_indi(outlier_indi, all_individuals)
		if cluster_id not in valid_clusters:
			continue
		t.write(cluster_id.split('_')[0] + '\t' + chrom_num + '\t' + variant_position + '\t' + ','.join(outlier_indi) + '\t' + ','.join(inlier_indi) + '\t' + distance_to_ss + '\t' + strand + '\t' + anno + '\t' + variant_allele + '\n')
	t.close()


cluster_info_file = sys.argv[1]
whole_blood_outlier_file = sys.argv[2]
xt_outlier_file = sys.argv[3]
watershed_posterior_file = sys.argv[4]
exon_file = sys.argv[5]
inlier_variant_position_file = sys.argv[6]
processed_data_output_dir = sys.argv[7]





cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)

step_1_output_file = processed_data_output_dir + 'background_step_1_filters.txt'
#run_step_1_filters(inlier_variant_position_file, cluster_to_strand_mapping, step_1_output_file)


step_2_output_file = processed_data_output_dir + 'background_step_2_filters.txt'
#run_step_2_filters(step_1_output_file, whole_blood_outlier_file, step_2_output_file, 'whole_blood_cluster_splicing_pvalue')
step_3_output_file = processed_data_output_dir + 'background_step_3_filters.txt'
#run_step_2_filters(step_2_output_file, xt_outlier_file, step_3_output_file, 'xt_cluster_splicing_pvalue')

# Add cluster info to file (ie add gene names and junction names)
step_4_output_file = processed_data_output_dir + 'background_step_4_filters.txt'
#run_step_4_add_cluster_info(step_3_output_file, cluster_info_file, step_4_output_file)

# Add Watershed INFO
step_5_output_file = processed_data_output_dir + 'background_step_5_filters.txt'
#run_step_5_add_watershed_info(step_4_output_file, watershed_posterior_file, step_5_output_file)

step_6_output_file = processed_data_output_dir + 'background_step_6_filters.txt'
#run_step_6_check_to_see_if_variant_has_high_watershed_in_any_individuals(step_5_output_file, watershed_posterior_file, step_6_output_file)


organized_output_file = processed_data_output_dir + 'background_variants.txt'
organized_with_junctions_output_file = processed_data_output_dir + 'background_variants_with_junc_positions.txt'
#organized_test_variants(step_6_output_file, organized_output_file, organized_with_junctions_output_file)

cluster_to_ss_mapping = extract_mapping_from_cluster_id_to_junctions(cluster_info_file)

cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)
# Create dictionary of annotated splice sites
# Where keys are "chr"$chrom_num"_"#splice_sites
annotated_splice_sites = get_dictionary_of_annotated_splice_sites(exon_file)




#create clusters_to_plot_file
clusters_to_plot_file = processed_data_output_dir + 'background_clusters_to_plot.txt'
make_clusters_to_plot_file(organized_output_file, cluster_to_ss_mapping, cluster_to_strand_mapping, annotated_splice_sites, whole_blood_outlier_file, clusters_to_plot_file)
