import numpy as np 
import os
import sys
import pdb

def merge_parallelized_results(output_root, suffix, total_jobs):
	# Open output (merged result) file handle
	t = open(output_root + 'merged_' + suffix, 'w')
	# Loop through parrallelized jobs
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + '_' + suffix
		# Open file for one job
		f = open(file_name)
		# To identify header
		head_count = 0
		# Stream file from one job
		for line in f:
			line = line.rstrip()
			# HEADER
			if head_count == 0:
				head_count = head_count + 1
				# Print header if this the first job
				if job_number == 0:
					t.write(line + '\n')
				continue
			# Standard line
			t.write(line + '\n')
		f.close()
		# Delete file from single job
		os.system ('rm ' + file_name)
	t.close()



output_root = sys.argv[1]
total_jobs = int(sys.argv[2])



merge_parallelized_results(output_root, "emperical_pvalue.txt", total_jobs)
merge_parallelized_results(output_root, "md.txt", total_jobs)