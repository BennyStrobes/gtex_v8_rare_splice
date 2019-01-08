import numpy as np 
import os
import sys
import pdb


#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file, output_file, missing_file, liftover_directory, direction):
	if direction == 'hg38_to_hg19':
		direction_word = 'hg38ToHg19'
	elif direction == 'hg19_to_hg38':
		direction_word = 'hg19ToHg38'
	else:
		print('assumption errror:' + direction + ' is not a valid argument')
		pdb.set_trace()
	stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + direction_word + '.over.chain.gz ' + output_file + ' ' + missing_file
	os.system(stringer)


input_bed_file = sys.argv[1]
output_bed_file = sys.argv[2]
liftover_directory = sys.argv[3]
liftover_direction = sys.argv[4]  # Either 'hg38_to_hg19' or 'hg19_to_hg38'

# Run liftover
temporary_missing_file = output_bed_file + '.missing'
run_liftover(input_bed_file, output_bed_file, temporary_missing_file, liftover_directory, liftover_direction)

# Remove temporary missing file
os.system('rm ' + temporary_missing_file)