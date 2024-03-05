#!/usr/bin/python

import os
import sys
import gzip
import numpy as np
import pandas as pd


##


# Args
path_input = sys.argv[1]
res = int(sys.argv[2])

# Paths
# path_input = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/functional_genomics/HiChip/filtered_loops/prova'
# res = i8

# Paths
filename = os.path.basename(path_input)
path_output = os.path.join(os.path.dirname(path_input), f'{filename}.tsv.gz')

# Open
filenames = sorted(os.listdir(path_input))
output = gzip.open(path_output, 'wt', encoding='utf-8')

# Stream
for i, x in enumerate(filenames):
    file_path = os.path.join(path_input, x)
    f = open(file_path, 'r', encoding='utf-8')
    if i == 0:
        output.writelines(f.readline())
    else:
        next(f)
    for line in f:
        output.write(line)
    f.close()

# Close
output.close()


##


# After concatenation, read, filter, subset and write clean file
dtypes = {
    'seqnames1':'str', 'start1':np.int64, 'end1':np.int64, 
    'seqnames2':'str', 'start2':np.int64, 'end2':np.int64, 
    'D':np.int64, 'counts':np.int16, 'mu':np.float16, 
    'sdev':np.float16, 'pvalue':np.float16, 'qvalue':np.float16
}
df = pd.read_csv(path_output, sep='\t', usecols=list(dtypes.keys()), dtype=dtypes)

# Filter significant loops and drop loops at outside the target distance range 
# and with too low (5) counts. Overwrite.
df.dropna().query('qvalue<=0.01').to_csv(path_output, sep='\t', index=False)


##