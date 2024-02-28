#!/usr/bin/python

import sys
import os
import numpy as np
import pandas as pd


# Args
path_input = sys.argv[1]
path_output = sys.argv[2]
filename = sys.argv[3]
res = int(sys.argv[4])

# Read
df = pd.read_csv(os.path.join(path_input, filename), sep='\t')

# Filter highly significative interactions 
thr = 3*res
counts = 10
df = df.query('qvalue<=0.01 and D>=@thr and D<=1000000 and counts>=5')
print(f'Filtered loops: {df.shape[0]}')

# Save
df.to_csv(os.path.join(path_output, filename), sep='\t', index=False)


##
