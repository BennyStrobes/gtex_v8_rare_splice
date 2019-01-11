import numpy as np 
import os
import sys
import pdb

# Get dictionary list of rare varants
def get_dictionary_list_of_rare_variants(variant_bed_file, position_to_ensambles):
	# Initialize dictionary
	dicti = {}
	# Stream variant bed file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent features
		individual_id = data[0]
		chrom_num = data[1]
		var_pos = data[2]
		maf = data[4]
		num_rv = data[5]
		ref = data[6]
		alt = data[7]
		#######################
		# NEED TO BE FIXED!!
		#######################
		if int(num_rv) == 2 or int(num_rv) == 0:
			continue
		if int(num_rv) > 2:
			pdb.set_trace()
		if int(num_rv) < 0:
			pdb.set_trace()
		if ref == alt:
			pdb.set_trace()

		# Skip variants not mapped to any gene
		if chrom_num + '_' + var_pos not in position_to_ensambles:
			continue
		gene_array = position_to_ensambles[chrom_num + '_' + var_pos]
		for gene in gene_array:
			variant_id = chrom_num + '_' + var_pos + '_' + individual_id + '_' + gene
			if variant_id not in dicti:
				dicti[variant_id] = 0
	f.close()
	return dicti

# Create mapping from genomic position to ensamble id
def get_mapping_from_position_to_ensamble_id(rare_variant_to_gene_file):
	# Initialize mapping
	mapping = {}
	# Stream input file 
	f = open(rare_variant_to_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relvent columns
		ensamble_id = data[1]
		chrom_num = data[2]
		chrom_pos = data[3]
		# Add to mapping
		variant_position_string = chrom_num + '_' + chrom_pos
		if variant_position_string not in mapping:
			mapping[variant_position_string] = {}
		mapping[variant_position_string][ensamble_id] = 0
	# Close stream
	f.close()
	return mapping

# Filtered lines in genomic annotation file to only include (variant, individual)
def filter_raw_genomic_annotation_file(raw_genomic_annotation_file, filtered_genomic_annotation_file, rare_variants):
	# Open input file handle
	f = open(raw_genomic_annotation_file)
	# Open output file handle
	t = open(filtered_genomic_annotation_file, 'w')
	# To skip header
	head_count = 0
	# Stream genomic annotation file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		# Extract relevent fields
		ensamble_id = data[0]
		indi_id = data[40].split('_')[0] + '-' + data[40].split('_')[1]
		chrom_num = data[41]
		var_pos = data[5]
		# Form variant id from relevent fields
		variant_id = chrom_num + '_' + var_pos + '_' + indi_id + '_' + ensamble_id
		# Throw out line if not in in rare_variants dictionary (ie not within 10KB from gene)
		if variant_id not in rare_variants:
			continue
		# print line to filtered output file
		t.write(line + '\n')
	# Close filehandles
	t.close()
	f.close()

# Compress lines from transcript level to gene level
def compress_annotations_from_transcript_level_to_gene_level(filtered_genomic_annotation_file, gene_level_genomic_annotation_file):
	# NOTES #
	# data[9] is 'ensembl_tok' and varies between transcripts
	# data[10] is 'parent_category' (created by Xin that varies) (TO BE REMOVED)
	# data[11] is 'LoF' which varies
	# data[12] is 'Feature_type' which varies
	# data[13] is 'Feature' which varies. This is transcript id ( TO BE REMOVED)
	# data[14] is 'BIOTYPE' which varies. TOOOO MANY. Probably remove
	# data[15] is 'ensembl_tss' which varies (created by Xin that varies) (TO BE REMOVED)

	# First pass through file to fill in dictionary mapping from (variant, individual, gene) to annotations
	dicti9 = {}
	dicti10 = {}
	dicti11 = {}
	dicti12 = {}
	dicti14 = {}
	dicti15 = {}
	# Used to skip header
	temp_dicti = {}
	head_count = 0
	f = open(filtered_genomic_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if len(data) != 43:
			print('length error')
			pdb.set_trace()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_id = data[0]
		var_pos = data[5]
		chrom_num = data[41]
		indi_id = data[40].split('_')[0] + '-' + data[40].split('_')[1]
		# Form variant id from relevent fields
		variant_id = chrom_num + '_' + var_pos + '_' + indi_id + '_' + gene_id
		# If never seen variant before add to dictionary
		dicti9[data[9]] = 1
		dicti10[data[10]] = 1
		dicti11[data[11]] = 1
		dicti12[data[12]] = 1
		dicti14[data[14]] = 1
		dicti15[data[15]] = 1

		

	f.close()
	#t.close()
	pdb.set_trace()
	return


raw_genomic_annotation_file = sys.argv[1]
variant_bed_file = sys.argv[2]
rare_variant_to_gene_file = sys.argv[3]
genomic_annotation_dir = sys.argv[4]

# Create mapping from genomic position to ensamble ids
#position_to_ensambles = get_mapping_from_position_to_ensamble_id(rare_variant_to_gene_file)


# Get dictionary list of rare varants
##############################################
# NEED TO FIX SKIP OF num_rv == 2 and num_rv == 0!!!!!
# Should also make sure alleleles in variant bed file match alleles in raw_genomic_annotation file
##############################################
#rare_variants = get_dictionary_list_of_rare_variants(variant_bed_file, position_to_ensambles)

# STEP 1
# Filtered lines in genomic annotation file to only include (variant, individual)
filtered_genomic_annotation_file = genomic_annotation_dir + 'filtered_raw_genomic_annotations.txt'
# filter_raw_genomic_annotation_file(raw_genomic_annotation_file, filtered_genomic_annotation_file, rare_variants)

# STEP 2
# Compress lines from transcript level to gene level
gene_level_genomic_annotation_file = genomic_annotation_dir + 'filtered_raw_gene_level_genomic_annotations.txt'
compress_annotations_from_transcript_level_to_gene_level(filtered_genomic_annotation_file, gene_level_genomic_annotation_file)