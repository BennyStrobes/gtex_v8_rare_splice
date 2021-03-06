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
		num_alt_alleles = data[5]
		ref = data[6]
		alt = data[7]

		# Skip variants not mapped to any gene
		if chrom_num + '_' + var_pos not in position_to_ensambles:
			continue
		gene_array = position_to_ensambles[chrom_num + '_' + var_pos]
		for gene in gene_array:
			variant_id = chrom_num + '_' + var_pos + '_' + individual_id + '_' + gene
			if variant_id not in dicti:
				dicti[variant_id] = ({ref:1,alt:1}, num_alt_alleles)
			else:
				if num_alt_alleles != dicti[variant_id][1]:
					print('assumption error ')
					pdb.set_trace()
				if ({ref:1,alt:1}) != dicti[variant_id][0]:
					print('assumption error')
					pdb.set_trace()
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
	used = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\tnum_rare_variants' + '\n')
			continue
		# Extract relevent fields
		ensamble_id = data[0]
		indi_id = data[40].split('_')[0] + '-' + data[40].split('_')[1]
		chrom_num = data[41]
		var_pos = data[5]
		ref = data[2]
		alt = data[3]
		rv_allele = data[4]
		# Form variant id from relevent fields
		variant_id = chrom_num + '_' + var_pos + '_' + indi_id + '_' + ensamble_id
		# Throw out line if not in in rare_variants dictionary (ie not within 10KB from gene)
		if variant_id not in rare_variants:
			continue
		# ERROR CHECKING
		alleles = rare_variants[variant_id][0]
		num_alt_alleles = rare_variants[variant_id][1]
		if num_alt_alleles == '1':
			if ref not in alleles:
				pdb.set_trace()
			if alt not in alleles:
				pdb.set_trace()
		elif num_alt_alleles == '0':
			if len(alleles) != 1:
				pdb.set_trace()
			if ref not in alleles:
				pdb.set_trace()
			if ref != rv_allele:
				pdb.set_trace()
		elif num_alt_alleles == '2':
			if len(alleles) != 1:
				pdb.set_trace()
			if alt not in alleles:
				pdb.set_trace()
			if alt != rv_allele:
				pdb.set_trace()
		else:
			print('what number of alleles')
			pdb.set_trace()
		# print line to filtered output file
		if num_alt_alleles == '0':
			num_alt_alleles_print = '2'
		else:
			num_alt_alleles_print = num_alt_alleles
		t.write(line + '\t' + num_alt_alleles_print + '\n')
	# Close filehandles
	t.close()
	f.close()


# Merge data to get gene level annotation
def merge_annotations_to_gene_level(data1, data2):
	# Error checking
	if np.array_equal(data1[:9],data2[:9]) == False:
		print('mismatch')
		pdb.set_trace()
	if np.array_equal(data1[11:], data2[11:]) == False:
		print('mismatch')
		pdb.set_trace()
	if len(data1) != len(data2):
		print('length mismatch')
		pdb.set_trace()
	new_column_9 = data1[9] + ',' + data2[9]
	new_column_10 = data1[10] + ',' + data2[10]
	return np.concatenate((data1[:9], [new_column_9], [new_column_10], data1[11:]))

# Compress lines from transcript level to gene level
def compress_annotations_from_transcript_level_to_gene_level(filtered_genomic_annotation_file, gene_level_genomic_annotation_file):
	# There are several features that vary for the same (individual, variant, gene). Thats because these features are at the transcript level.
	# They are:
	# data[9]: 'ensembl_tok' and varies between transcripts (take union)
	# data[10]:'parent_category' (created by Xin from ensamble_tok) (TO BE REMOVED)
	# data[11]:'LoF' (take union)
	# data[12]: 'Feature_type' (TO BE REMOVED)
	# data[13]: 'Feature' which varies. This is transcript id ( TO BE REMOVED)
	# data[14]: BIOTYPE'. TOOOO MANY ( TO BE REMOVED)
	# data[15]: ensembl_tss' (created by Xin from ensamble_tok) (TO BE REMOVED)

	# Columns to be kept
	valid_columns = range(10) + [11] + range(16,44)

	# First pass through file to fill in dictionary mapping from (variant, individual, gene) to annotations
	# Initialize dictionary 
	dicti = {}
	# Used to skip header
	head_count = 0
	# Stream file
	f = open(filtered_genomic_annotation_file)
	counter = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split('\t'))
		# Error checking
		if len(data) != 44:
			print('length error')
			pdb.set_trace()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			header = data[valid_columns]
			continue
		counter = counter + 1
		if counter%1000000 == 0:
			print(counter)
		# Extract relevent fields
		gene_id = data[0]
		var_pos = data[5]
		chrom_num = data[41]
		indi_id = data[40].split('_')[0] + '-' + data[40].split('_')[1]
		# Form variant id from relevent fields
		variant_id = chrom_num + '_' + var_pos + '_' + indi_id + '_' + gene_id
		# If never seen variant before add to dictionary
		if variant_id not in dicti:
			# Remove columns 10, 12, 13, 14, 15
			filtered_data = data[valid_columns]
			# Add to dictionary
			dicti[variant_id] = filtered_data
		# We have seen this variant before in another transcript
		else:
			# Now we have to merge annotations
			# Remove columns 10, 12, 13, 14, 15
			filtered_data = data[valid_columns]
			other_filtered_data = dicti[variant_id]
			# Merge data to get gene level annotation
			merged_data = merge_annotations_to_gene_level(filtered_data, other_filtered_data)
			dicti[variant_id] = merged_data
	f.close()
	# Open output file
	t = open(gene_level_genomic_annotation_file, 'w')
	# Print header to output file
	t.write('\t'.join(header) + '\n')
	# Loop through (variant, gene, individuals)
	for variant in dicti.keys():
		data = dicti[variant]
		t.write('\t'.join(data[:9]) + '\t')
		# Get unique elements of data[9] and data[10]
		t.write(','.join(np.unique(data[9].split(','))) + '\t')
		t.write(','.join(np.unique(data[10].split(','))) + '\t')
		t.write('\t'.join(data[11:]) + '\n')
	# close output file
	t.close()
	return

# STEP 3
# Make genomic annotations real valued
# Many of the genomic annotations are categorical. Convert these categorical variables to sets of binary variables.
def make_real_valued_genomic_annotations(input_file, output_file):
	# NOTES:
	# 0: ensamble ID (always starts with 'ENSG') (No Nans)
	# 1: Chrom_numeric  (No Nans)
	# 2: Ref (No Nans)
	# 3: Alt (No Nans)
	# 4: allele (No Nans)
	# 5: position (No Nans)
	# 6: gene type (No nans) (Too many clases; TO BE REMOVED)
	# 7: distTSS (No Nans)
	# 8: distTES (No Nans)
	# 9: ensamble_tok: Categorical. Currently seperated by commas. Definitely some NAs
	# 10: LoF: Categorical: Currently seperated by commas. Some NA and <undefined>
	# 11. VARIANT_CLASS: Currently SNV or undefined (Too few classes. all snvs; TO BE REMOVED)
	# 12. af_1kg: (Some NaN and numeric) (TO BE REMOVED because data[33] is better)
	# 13. phylop: Some NaN and numeric (impute: 0.0) (not in cadd, but used approximate mean)
	# 14. gc: Some NaN and numeric (impute: .418)
	# 15. cpg: Some NaN and numeric (impute: .024)
	# 16. siftcat: Some 'NA' and '<undefined>'. Otherwise takes on 'deleterious' and 'tolerated'
	# 17. siftval: NaN and numeric (impute: 0)
	# 18. polyphencat: Categorical: Some '<undefined>' and 'NA'
	# 19. polyphenval: NaN and numeric (impute: 0)
	# 20. bstatistic: NaN and numeric (800.261)
	# 21. priphcons: NaN and numeric (.115)
	# 22. mamphcons: NaN and numeric (.079)
	# 23. verphcons: NaN and numeric (.094)
	# 24. priphylop: NaN and numeric (-.033)
	# 25. mamphylop: NaN and numeric (-.038)
	# 26. verphylop: NaN and numeric (.017)
	# 27. gerpn: NaN and numeric (1.909)
	# 28. gerps: NaN and numeric (-.2)
	# 29: phred: NaN and numeric (0.0. not in cad but ok)
	# 30: indiv_id: (TO BE REMOVED)
	# 31: GT: Is only '1|0' or '0|1' (TO BE REMOVED)
	# 32: genotype: (To Be removed)
	# 33: af (No Nans)
	# 34: RACE: As we already filtered on european is only 'White' (TO BE REMOVED)
	# 35: otherId: Gtex Id
	# 36: chr: chrom num (No NAns)
	# 37: af_gnomad: Some NAs. (data[33] is better here) (TO BE REMOVED)
	###############################

	#########################
	# Create Imputation dictionary
	impute = {}
	impute[13] = 0.0
	impute[14] = .418
	impute[15] = .024
	impute[17] = 0
	impute[19] = 0
	impute[20] = 800.261
	impute[21] = .115
	impute[22] = .079
	impute[23] = .094
	impute[24] = -.033
	impute[25] = -.038
	impute[26] = .017
	impute[27] = 1.909
	impute[28] = -.2
	impute[29] = 0.0

	#########################
	# Create dictionary containing columns to be removed
	to_be_removed = {}
	to_be_removed[0] = 1
	to_be_removed[1] = 1
	to_be_removed[2] = 1
	to_be_removed[3] = 1
	to_be_removed[4] = 1
	to_be_removed[5] = 1
	to_be_removed[6] = 1
	to_be_removed[11] = 1
	to_be_removed[12] = 1
	to_be_removed[30] = 1
	to_be_removed[31] = 1
	to_be_removed[32] = 1
	to_be_removed[34] = 1
	to_be_removed[35] = 1
	to_be_removed[36] = 1
	to_be_removed[37] = 1

	#########################
	# Create dictionary containing that need no change
	no_change = {}
	no_change[7] = 1
	no_change[8] = 1
	no_change[33] = 1
	no_change[38] = 1
	no_change[39] = 1
	no_change[40] = 1
	no_change[41] = 1
	no_change[42] = 1
	no_change[43] = 1


	#########################
	# Create dictionary containing columns that need to be converted from categorical into binaries
	categoricals = {}
	categoricals[9] = []
	categoricals[10] = []
	categoricals[16] = []
	categoricals[18] = []

	for index in range(44,71):
		categoricals[index] = []


	# First loop through file one time to fill in cateogricals dictionary
	# For header
	head_count = 0
	# Stream input file
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Header 
		if head_count == 0:
			head_count = head_count + 1
			header = data
			continue
		for i, ele in enumerate(data):
			if i in categoricals:
				ele_array = ele.split(',')
				for small_ele in ele_array:
					if small_ele == '<undefined>' or small_ele == 'NA':
						continue
					categoricals[i].append(small_ele)
	f.close()
	# Take only uniques
	for key in categoricals.keys():
		categoricals[key] = np.unique(categoricals[key])
	head_count = 0
	print('start')
	# Open output file handle
	t = open(output_file, 'w')
	# Take another pass through
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('SubjectID\tGeneName\tvariant_chromosome\tvariant_position\tvariant_id')
			for i,ele in enumerate(data):
				if i in to_be_removed:
					continue
				elif i in no_change:
					t.write('\t' + ele)
				elif i in impute:
					t.write('\t' + ele)
				elif i in categoricals:
					for namer in categoricals[i]:
						t.write('\t' + ele + '_' + namer)
			t.write('\n')
			continue
		# Id field
		rv_allele = data[4]
		alt = data[3]
		ref = data[2]
		num_rv = data[38]
		if rv_allele == alt:
			major_allele = ref 
		else:
			major_allele = alt
		variant_id = data[36] + '_' + data[5] + '_' + major_allele + '_' + rv_allele + '_' + num_rv
		t.write(data[35] + '\t' + data[0] + '\t' + data[36] + '\t' + data[5] + '\t' + variant_id)
		for i, ele in enumerate(data):
			if i in to_be_removed:
				continue
			elif i in no_change:
				t.write('\t' + ele)
			elif i in impute:
				if ele.startswith('Na') or ele.startswith('NA'):
					t.write('\t' + str(impute[i]))
				else:
					t.write('\t' + ele)
			elif i in categoricals:
				initial_array = np.zeros(len(categoricals[i]))
				for itera,namer in enumerate(categoricals[i]):
					if namer == ele:
						initial_array[itera] = 1.0
				t.write('\t' + '\t'.join(initial_array.astype(str)))
			else:
				print('assumption error')
				pdb.set_trace()
		t.write('\n')
	f.close()
	t.close()

def compress_feature_vec_of_multiple_variants_for_the_same_gene(vec1, vec2):
	if len(vec1) != len(vec2):
		print('length assumption error')
		pdb.set_trace()
	part1 = np.minimum(vec1[:2],vec2[:2])
	part2 = np.maximum(vec1[2:47], vec2[2:47])
	part3 = np.asarray([np.minimum(vec1[47],vec2[47])])
	part4 = np.asarray([vec1[48] + vec2[48]])
	part5 = np.maximum(vec1[49:], vec2[49:])
	return np.concatenate((part1,part2,part3, part4,part5))


# STEP 4
# Compress multiple variants onto same gene
def compress_multiple_variants_onto_same_gene(input_file, output_file, variants_mapped_to_gene_file):
	# Most annotations we take maximum over.
	# The following we take the minimum over:
	# 1. distTSS
	# 2. distTES
	# 3. af
	# Make dictionary where keys are (individual, gene) pairs and values are vectors corresponding to their features
	mapping = {}
	variant_mapping = {}
	# Stream input file handle
	# Used to get header
	head_count = 0
	f = open(input_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header_features = data[5:]
			continue
		# Extract relevent fields from line
		indi_id = data[0].split('_')[0] + '-' + data[0].split('_')[1]
		ensamble_id = data[1]
		pair_name = indi_id + '_' + ensamble_id
		variant_id = data[4]
		feature_vec = np.asarray(data[5:]).astype(float)
		# Make distances absolute value
		feature_vec[0] = np.abs(feature_vec[0])
		feature_vec[1] = np.abs(feature_vec[1])
		# Add pair
		if pair_name not in mapping:
			mapping[pair_name] = feature_vec
		else:
			mapping[pair_name] = compress_feature_vec_of_multiple_variants_for_the_same_gene(feature_vec, mapping[pair_name])
		if pair_name not in variant_mapping:
			variant_mapping[pair_name] = []
		variant_mapping[pair_name].append(variant_id)
	f.close()
	# Open output file
	t = open(output_file, 'w')
	# Print header
	t.write('SubjectID\tGeneName\t' + '\t'.join(header_features) + '\n')
	for pair_name in mapping.keys():
		indi_id = pair_name.split('_')[0]
		ensamble_id = pair_name.split('_')[1]
		t.write(indi_id + '\t' + ensamble_id + '\t')
		feature_vec = mapping[pair_name].astype(str)
		t.write('\t'.join(feature_vec) + '\n')
	t.close()
	t = open(variants_mapped_to_gene_file, 'w')
	t.write('SubjectID\tGeneName\tvariants\n')
	for pair_name in variant_mapping.keys():
		indi_id = pair_name.split('_')[0]
		ensamble_id = pair_name.split('_')[1]
		t.write(indi_id + '\t' + ensamble_id + '\t')
		variant_vec = np.unique(variant_mapping[pair_name]).astype(str)
		t.write(','.join(variant_vec) + '\n')
	t.close()

def get_cadd_score(chrom_num, pos, major_allele, rv_allele, cadd_score_tabix):
	# Binary variable on whether our Rare variant is scored by cadd (initialized to False)
	rv_scored_by_cadd = False
	# Initialize cadd score to imputed value of score
	score = '1.0'
	# Looop through all variants in this position in CADD
	for poss_var in cadd_score_tabix.fetch(chrom_num,pos-1,pos):
		cadd_data = poss_var.rstrip().split('\t')
		cadd_ref = cadd_data[2]
		cadd_alt = cadd_data[3]
		phred_cadd_score = cadd_data[5]
		if cadd_ref == major_allele and cadd_alt == rv_allele:
			if rv_scored_by_cadd == True:
				print('repeat')
				pdb.set_trace()
			rv_scored_by_cadd = True
			score = phred_cadd_score
	if rv_scored_by_cadd == False:
		print('miss')
	return score

# Step 2a
# Add annotations to genomic file
def add_annotations(gene_level_genomic_annotation_file, gene_level_genomic_annotation_add_anno_file, cadd_file, cadd_anno_file):
	# Open Tabix files for genomic annotations
	cadd_score_tabix = pysam.Tabixfile(cadd_file,'r') 
	cadd_anno_tabix = pysam.Tabixfile(cadd_anno_file, 'r')
	# Stream input file
	f = open(gene_level_genomic_annotation_file)
	# Used to skip header
	head_count = 0
	count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		chrom_num = data[36]
		chrom_num = data[1]
		pos = int(data[5])
		ensamble_id = data[0]
		ref = data[2]
		alt = data[3] 
		rv_allele = data[4]
		if rv_allele == alt:
			major_allele = ref 
		else:
			major_allele = alt
			print('potential miss')

		cadd_score = get_cadd_score(chrom_num, pos, major_allele, rv_allele, cadd_score_tabix)
		count = count + 1
		print(count)

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

# Create mapping from cluster id to strand
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

def get_gene_to_clusters_mapping(cluster_info_file):
	# Used to skip header
	head_count = 0
	# Initialize mapping 
	mapping = {}
	# Stream cluster info file
	f = open(cluster_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_array = data[2].split(',')
		cluster_id = data[0]
		# For each gene associated to this cluster
		for gene in gene_array:
			# Add gene to mapping if not already there
			if gene not in mapping:
				mapping[gene] = []
			# Add cluster to gene's array
			mapping[gene].append(cluster_id)
	f.close()
	# Mapping from gene to unique set of clusters
	for gene_id in mapping.keys():
		mapping[gene_id] = np.unique(mapping[gene_id])
	return mapping

def get_cluster_level_ss_feature_vector(cluster_mapping, var_pos, chrom_num, strand):
	donor_variant = 0
	acceptor_variant = 0
	ppt_variant = 0
	donor_window_variant = 0
	acceptor_window_variant = 0
	for start_splice_site in cluster_mapping['start_splice_sites']:
		distance = start_splice_site - var_pos
		# Check if its a donor ss
		if strand == '+': # Donor ss
			if distance >= -6 and distance <= 0:
				donor_variant = 1
			if distance >= -10 and distance < -6:
				donor_window_variant = 1
			if distance > 0 and distance <= 10:
				donor_window_variant = 1
		elif strand == '-': # Acceptor ss
			if distance >= -2 and distance <= 0:
				acceptor_variant = 1
			if distance <= -5 and distance > -35:
				ppt_variant = 1
			if distance >= -4 and distance < -2:
				acceptor_window_variant = 1
			if distance > 0 and distance <= 10:
				acceptor_window_variant = 1
	for end_splice_site in cluster_mapping['end_splice_sites']:
		distance = var_pos - end_splice_site
		if strand == '-': # Donor ss
			if distance >= -6 and distance <= 0:
				donor_variant = 1
			if distance >= -10 and distance < -6:
				donor_window_variant = 1
			if distance > 0 and distance <= 10:
				donor_window_variant = 1
		elif strand == '+': # Acceptor ss
			if distance >= -2 and distance <= 0:
				acceptor_variant = 1
			if distance <= -5 and distance > -35:
				ppt_variant = 1
			if distance >= -4 and distance < -2:
				acceptor_window_variant = 1
			if distance > 0 and distance <= 10:
				acceptor_window_variant = 1
	return np.asarray([donor_variant, acceptor_variant, ppt_variant, donor_window_variant, acceptor_window_variant])

# Add observed splice site information
def add_observed_ss_genomic_annotations(input_file, output_file,exon_file, cluster_info_file):
	# Create mapping from cluster id to junctions
	cluster_mapping = extract_mapping_from_cluster_id_to_junctions(cluster_info_file)
	# Create mapping from cluster id to strand
	cluster_to_strand_mapping = get_cluster_to_strand_mapping(cluster_info_file, exon_file)
	# Create mapping from gene to clusters
	gene_to_clusters_mapping = get_gene_to_clusters_mapping(cluster_info_file)
	# Open input and output filehandles
	f = open(input_file)
	t = open(output_file, 'w')
	# Used to skip header
	count = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\tdonor_ss\tacceptor_ss\tppt_region\tdonor_ss_window\tacceptor_ss_window\n')
			continue
		# Extract relevent fields
		chrom_num = data[36]
		ensamble_id = data[0] 
		position = int(data[5])
		# Use mapping to extract list of clusters associated with this gene
		if ensamble_id in gene_to_clusters_mapping:  # Jxn data for gene
			clusters = gene_to_clusters_mapping[ensamble_id]
			# Initialize ss annotation vec for this gene
			ss_gene_level_vec = np.asarray([0,0,0,0,0])
			for cluster in clusters: 
				# Extract splice site feature vector for this cluster
				ss_vec = get_cluster_level_ss_feature_vector(cluster_mapping[cluster], position, chrom_num, cluster_to_strand_mapping[cluster])
				# Take maximum across clusters (Does the observed ss occur in ANY cluster for that gene)
				ss_gene_level_vec = np.maximum(ss_vec, ss_gene_level_vec)
		else:  # No jxn data for gene
			ss_gene_level_vec = np.asarray([0,0,0,0,0])
		t.write(line + '\t' + '\t'.join(ss_gene_level_vec.astype(str)) + '\n')
	f.close()
	t.close()

def clean_up_gene_level_real_valued_genomic_annotations(input_file, output_file):
	f = open(input_file)
	t = open(output_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(data[0] + '\t' + data[1] + '\t' + '\t'.join(data[5:]) + '\n')
			continue
		indi_id = data[0].split('_')[0] + '-' + data[0].split('_')[1]
		d1 = str(abs(float(data[5])))
		d2 = str(abs(float(data[6])))
		t.write(indi_id + '\t' + data[1] + ':' + data[4] + '\t' + d1 + '\t' + d2 + '\t' + '\t'.join(data[7:]) + '\n')
	f.close()
	t.close()

def convert_to_broader_chrom_hmm_state(arr, mapping):
	new_arr = []
	for ele in arr:
		words = []
		small_arr = ele.split(',')
		for small_ele in small_arr:
			words.append(mapping[small_ele])
		new_arr.append(','.join(words))
	return new_arr

# STEP 3b
# Add tissue-specic chromHMM elements
def add_tissue_specific_chromhmm_annotations(input_file, output_file, chrom_hmm_file):
	# Create mapping from chromHMM state to broader state
	mapping = {}
	for chrom_hmm_state in range(1,16):
		if chrom_hmm_state == 6 or chrom_hmm_state == 7 or chrom_hmm_state == 12:
			mapping[str(chrom_hmm_state)] = 'enh'
		elif chrom_hmm_state == 1 or chrom_hmm_state == 2 or chrom_hmm_state == 10:
			mapping[str(chrom_hmm_state)] = 'prom'
		elif chrom_hmm_state == 11:
			mapping[str(chrom_hmm_state)] = 'prom,enh'
		else:
			mapping[str(chrom_hmm_state)] = 'NA'
	mapping['.'] = 'NA'
	# Create dictionary mapping each variant to its broader chromHMM state
	f = open(chrom_hmm_file)
	head_count = 0
	counter = 0
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			header = []
			for ele in data[2:]:
				header.append('chrom_hmm_' + ele)
			continue
		variant_name = data[0] + '_' + data[1]
		dicti[variant_name] = convert_to_broader_chrom_hmm_state(data[2:], mapping)
		counter = counter + 1
	f.close()
	print('start')
	null_vec = ['NA']*len(header)
	f = open(input_file)
	t = open(output_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + '\t'.join(header) + '\n')
			continue
		t.write(line + '\t')
		variant_name = data[36] + '_' + data[5]
		if variant_name in dicti:
			t.write('\t'.join(dicti[variant_name]) + '\n')
		else:
			t.write('\t'.join(null_vec) + '\n')
	f.close()
	t.close()


raw_genomic_annotation_file = sys.argv[1]
variant_bed_file = sys.argv[2]
rare_variant_to_gene_file = sys.argv[3]
genomic_annotation_dir = sys.argv[4]
exon_file = sys.argv[5]
cluster_info_file = sys.argv[6]
chrom_hmm_file = sys.argv[7]


# Create mapping from genomic position to ensamble ids
position_to_ensambles = get_mapping_from_position_to_ensamble_id(rare_variant_to_gene_file)


# Get dictionary list of rare varants
##############################################
rare_variants = get_dictionary_list_of_rare_variants(variant_bed_file, position_to_ensambles)

# STEP 1
# Filtered lines in genomic annotation file to only include (variant, individual)
# NOTE: Every variant-gene-individual from rare_variants is present in the file
filtered_genomic_annotation_file = genomic_annotation_dir + 'filtered_raw_genomic_annotations.txt'
filter_raw_genomic_annotation_file(raw_genomic_annotation_file, filtered_genomic_annotation_file, rare_variants)

# STEP 2
# Compress lines from transcript level to gene level
gene_level_genomic_annotation_file = genomic_annotation_dir + 'filtered_raw_gene_level_genomic_annotations.txt'
compress_annotations_from_transcript_level_to_gene_level(filtered_genomic_annotation_file, gene_level_genomic_annotation_file)

# STEP 3a
# Add observed splice site information
gene_level_genomic_annotation_and_ss_file = genomic_annotation_dir + 'filtered_raw_gene_level_and_ss_genomic_annotations.txt'
add_observed_ss_genomic_annotations(gene_level_genomic_annotation_file, gene_level_genomic_annotation_and_ss_file,exon_file, cluster_info_file)

# STEP 3b
# Add tissue-specic chromHMM elements
gene_level_genomic_annotation_and_ss_and_chromhmm_file = genomic_annotation_dir + 'filtered_raw_gene_level_and_ss_and_chromHMM_genomic_annotations.txt'
add_tissue_specific_chromhmm_annotations(gene_level_genomic_annotation_and_ss_file, gene_level_genomic_annotation_and_ss_and_chromhmm_file, chrom_hmm_file)


# STEP 4a
# Make genomic annotations real valued
# Many of the genomic annotations are categorical. Convert these categorical variables to sets of binary variables.
gene_level_real_valued_genomic_annotation_file = genomic_annotation_dir + 'filtered_real_valued_gene_level_genomic_annotations.txt'
make_real_valued_genomic_annotations(gene_level_genomic_annotation_and_ss_and_chromhmm_file, gene_level_real_valued_genomic_annotation_file)

# STEP 4b
# Make file that has line for each variant in the same format as step 5.
gene_level_real_valued_clean_genomic_annotation_file = genomic_annotation_dir + 'filtered_real_valued_gene_level_clean_genomic_annotations.txt'
clean_up_gene_level_real_valued_genomic_annotations(gene_level_real_valued_genomic_annotation_file, gene_level_real_valued_clean_genomic_annotation_file)

print('mid')
# STEP 5
# Compress multiple variants onto same gene
gene_level_real_valued_variant_compressed_genomic_annotation_file = genomic_annotation_dir + 'filtered_real_valued_gene_level_variant_compressed_genomic_annotations.txt'
variants_mapped_to_gene_file = genomic_annotation_dir + 'variants_compressed_onto_genes.txt'
compress_multiple_variants_onto_same_gene(gene_level_real_valued_genomic_annotation_file, gene_level_real_valued_variant_compressed_genomic_annotation_file, variants_mapped_to_gene_file)
