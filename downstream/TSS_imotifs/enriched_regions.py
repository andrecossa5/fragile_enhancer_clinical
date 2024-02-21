"""
Find enriched regions into a FASTA file.
"""

import os
import pandas as pd
import numpy as np


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/TSS_motifs'
path_data = os.path.join('data', 'TSS_motifs')
path_main = os.path.join('results', 'TSS_motifs')


# Read FASTA
plus_seqs = []
f = open(os.path.join(path_data, 'FASTAfromBLISSends_FASTAminus.txt'), 'r')
for s in f.readlines():
    if not s.startswith('>'):
        plus_seqs.append(list(s.upper().strip('\n')))
f.close()

df = pd.DataFrame(plus_seqs)
top_represented = []
for i in range(df.shape[1]):
    top_represented.append(df.iloc[:,i].value_counts().index[0])


''.join(top_represented)

# 'CCCCCCCCCGGGTGGGCCCCC'
# 'GGGGGCCCACCCGGGCGGGGG'
