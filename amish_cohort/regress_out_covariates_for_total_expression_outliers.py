import numpy as np 
import os
import sys
import pdb
from sklearn import linear_model







####################
# command line args
####################
input_file = sys.argv[1]
output_file = sys.argv[2]

num_expression_pcs = 15

# Load in expression data
expr_full = np.loadtxt(input_file,dtype=str,delimiter='\t')
samples = expr_full[0,1:]
genes = expr_full[1:,0]
expr = np.transpose(expr_full[1:,1:].astype(float))

# RUN PCA
uuu, sss, vh = np.linalg.svd(np.transpose(expr))
expr_pc_loadings = np.transpose(vh)[:,:num_expression_pcs]

# Fit linear model
model = linear_model.LinearRegression(fit_intercept=True)
modelfit = model.fit(expr_pc_loadings,expr)
pred = modelfit.predict(expr_pc_loadings)

# Compute residual
resid = expr - pred


# Print output
f = open(input_file)
t = open(output_file,'w')

head_count = 0
counter = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	resid_expr = resid[:, counter]
	t.write(data[0] + '\t' + '\t'.join(resid_expr.astype(str)) + '\n')
	counter = counter + 1
f.close()
t.close()