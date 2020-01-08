import numpy as np 
import os
import sys
import pdb
import scipy.stats





####################
# command line args
####################
zscore_file = sys.argv[1]
pvalue_file = sys.argv[2]


head_count = 0
f = open(zscore_file)
t = open(pvalue_file, 'w')

for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	ensamble_id = data[0]
	zscores = np.asarray(data[1:]).astype(float)
	t.write(ensamble_id)
	for zscore in zscores:
		pvalue = scipy.stats.norm.sf(abs(zscore))*2
		t.write('\t' + str(pvalue))
	t.write('\n')
f.close()
t.close()
