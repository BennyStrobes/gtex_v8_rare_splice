import numpy as np 
import os
import sys
import pdb
import gzip








gtex_watershed_file = sys.argv[1]
gtex_variant_position_file = sys.argv[2]

t = open(gtex_variant_position_file, 'w')
f = gzip.open(gtex_watershed_file)

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		continue
	var = data[0].split(':')[2]
	chrom = var.split('_')[0]
	pos = var.split('_')[1]
	t.write(chrom + '\t' + pos + '\n')
f.close()
t.close()
