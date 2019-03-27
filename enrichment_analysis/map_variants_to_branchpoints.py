import numpy as np 
import os
import sys
import pdb


def make_branchpoint_chromosome(chrom_num_int, branchpoint_bed_file, distance_window):
	chrom = [0]*259250621
	chrom_string = 'chr' + str(chrom_num_int)
	f = open(branchpoint_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		line_chrom = data[0]
		if line_chrom != chrom_string:
			continue
		branchpoint = int(data[1])
		positions = range(branchpoint-distance_window, branchpoint+distance_window+1)
		for position in positions:
			chrom[position] = 1
	f.close()
	return chrom

variant_bed_file = sys.argv[1]  # Input file containing rare variants mapped to genes
branchpoint_bed_file = sys.argv[2]  #  Input file containing beds of branchpoint windows
variant_bed_branchpoint_mapped = sys.argv[3]  # Output file where variant file contains a column indicating whether variant overlaps branchpoint window
distance_window = int(sys.argv[4])

# Open output file 
t = open(variant_bed_branchpoint_mapped, 'w')

# Perform analysis chromosome by chromosome
for chrom_num_int in range(1,23):
	# Make array of length of chromosome. Zero if doesn't overlap brnachpoint window. 1 if it does
	chromosome = make_branchpoint_chromosome(chrom_num_int, branchpoint_bed_file, distance_window)

	chrom_string = 'chr' + str(chrom_num_int)
	# Stream variant bed file
	f = open(variant_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		line_chrom = data[0]
		if line_chrom != chrom_string:
			continue
		position = int(data[1])
		binary_fall_in_branchpoint = chromosome[position]
		t.write(line + '\t' + str(binary_fall_in_branchpoint) + '\n')
	f.close()
	t.flush()
t.close()


