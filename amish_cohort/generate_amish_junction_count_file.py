import numpy as np 
import os
import sys
import pdb
import gzip




def extract_gtex_junctions(gtex_lymphocyte_jxn_count_file):
	f = open(gtex_lymphocyte_jxn_count_file)
	head_count = 0
	dicti = {}
	array = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[0])
		junction_name = data[0].split(':')[0] + ':' + data[0].split(':')[1] + ':' + data[0].split(':')[2]
		dicti[junction_name] = {}
		# Error checking
		start = int(data[0].split(':')[1])
		end = int(data[0].split(':')[2])
		if end < start:
			print('assumption erroror')
			pdb.set_trace()
	return np.asarray(array), dicti

def fill_in_amish_sample_junction_counts(sample_name, junction_file, junction_dictionary):
	# Keep track of how many times junction occurs
	counts = {}
	f = gzip.open(junction_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		chrom_num = 'chr' + data[0]
		start_pos = str(int(data[1]))
		end_pos = str(int(data[2])+1)
		junction = chrom_num + ':' + start_pos + ':' + end_pos
		jxn_counts = int(data[4])

		# Simple error checking
		if data[3] != '.':
			print('assumption errro')
			pdb.set_trace()
		if float(end_pos) < float(start_pos):
			print('assumption error')
		if junction not in counts:
			counts[junction] = 0
		counts[junction] = counts[junction] + 1

		# Check to see if junction in junction file
		if junction not in junction_dictionary:
			continue
		if sample_name not in junction_dictionary[junction]:
			junction_dictionary[junction][sample_name] = 0
		junction_dictionary[junction][sample_name] = junction_dictionary[junction][sample_name] + jxn_counts
	f.close()

	# Simple error checking
	for junction in counts.keys():
		if counts[junction] > 2:
			print('assumption eorroro')
			pdb.set_trace()
	return junction_dictionary

def extract_amish_junction_counts(individual_to_junction_file_list, junction_dictionary):
	sample_names = []
	# Loop through file that has sample name along with junction file
	f = open(individual_to_junction_file_list)
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		sample_name = data[0]
		sample_names.append(sample_name)

		junction_file = data[1]
		#print(sample_name)
		# Fill in junction counts
		junction_dictionary = fill_in_amish_sample_junction_counts(sample_name, junction_file, junction_dictionary)
	return np.asarray(sample_names), junction_dictionary










#####################
# Command line args
#####################
gtex_lymphocyte_jxn_count_file = sys.argv[1] # File telling us what junctions to look for
individual_to_junction_file_list = sys.argv[2]
amish_junction_count_file = sys.argv[3] # output file


# Extract junctions used in gtex lymphocytes (these will be the only junctions we extract from the amish cohort)
ordered_junctions, junction_dictionary = extract_gtex_junctions(gtex_lymphocyte_jxn_count_file)

# Extract counts of the above junctions based on amish samples
ordered_samples, junction_dictionary = extract_amish_junction_counts(individual_to_junction_file_list, junction_dictionary)



# Print to output file 
t = open(amish_junction_count_file,'w')
# Print header
t.write('chrom' + '\t' + '\t'.join(ordered_samples) + '\n')

for full_junction_name in ordered_junctions:
	short_junction_name = full_junction_name.split(':')[0] + ':' + full_junction_name.split(':')[1] + ':' + full_junction_name.split(':')[2]
	t.write(full_junction_name)
	for i, sample_name in enumerate(ordered_samples):
		if sample_name in junction_dictionary[short_junction_name]:
			t.write('\t' + str(junction_dictionary[short_junction_name][sample_name]))
		else:
			t.write('\t0')
	t.write('\n')

t.close()
