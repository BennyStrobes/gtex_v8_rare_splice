import numpy as np 
import os
import sys
import pdb










#######################
# command line args
#######################
input_file = sys.argv[1]
output_file = sys.argv[2]



f = open(input_file)
t = open(output_file, 'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	ensamble_id = data[0]
	expression = np.asarray(data[1:]).astype(float)
	standardized_expression = (expression - np.mean(expression))/np.std(expression)
	t.write(ensamble_id + '\t' + '\t'.join(standardized_expression.astype(str)) + '\n')
f.close()
t.close()