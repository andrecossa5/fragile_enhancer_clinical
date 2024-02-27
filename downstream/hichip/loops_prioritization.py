"""
Visualize and quantify overlap between loops at different resolution.
"""

import os
import numpy as np
import pandas as pd
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Utils
def filter_enhancers(enh, df):
    """
    Test if as set of enhancers coincide with some loop interaction pair.
    """
    L = []
    for i in range(enh.shape[0]):
        d = enh.iloc[i].to_dict()
        enh_chr = d['chrI']
        enh_pos = d['startI']
        cluster = d['cluster']
        testI = (enh_pos >= df['startI']) & (enh_pos <= df['endI']) & (enh_chr == df['chrI'])
        testJ = (enh_pos >= df['startJ']) & (enh_pos <= df['endJ']) & (enh_chr == df['chrJ'])
        if testI.any():
            isin = 0
            idx = np.where(testI)[0]
        elif testJ.any():
            isin = 1
            idx = np.where(testJ)[0]
        start_i, end_i = df.iloc[idx,1:3].values[0]
        start_j, end_j = df.iloc[idx,4:6].values[0]
        l = [enh_chr, start_i, end_i, enh_chr, start_j, end_j, enh_pos, cluster, isin]
        L.append(l)
    
    # DataFrame
    df_coords = pd.DataFrame(L, columns=['chr', 'startI', 'endI', 'chr', 'startJ', 'endJ', 'enhancer_summit', 'cluster', 'isin_bin'])
    
    return df_coords
        

##


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical'
path_results = os.path.join(path_main, 'results', 'Hi_Chip')
path_enhancers = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/Hi_Chip/CtIP_enhancers.txt'
res = 8


##


# Viz per sample
scr = pd.read_csv(os.path.join(path_results, 'loops', 'SCR', f'SCR_loops_{res}.txt.gz'), sep='\t')
kd = pd.read_csv(os.path.join(path_results, 'loops', 'KD', f'KD_loops_{res}.txt.gz'), sep='\t')

# Filter highly significative interactions at proper distance (i.e., D >= 3*res*resolution and below 1MB)
thr = 3*1000*res
scr = scr.query('qvalue<=0.01 and D>=@thr and D<=1000000')
kd = kd.query('qvalue<=0.01 and D>=@thr and D<=1000000')

# Intersections and rankings based on loop strenghts (estimated mu)

# KD
kd_specific = kd.merge(scr, how='left', on=['chrI', 'startI', 'endI', 'chrJ', 'startJ', 'endJ'], suffixes=['_KD', '_SCR'])
test = kd_specific.loc[:,kd_specific.columns.str.contains('SCR')].isna().sum(axis=1) != 0
kd_specific = kd_specific.loc[test,:].iloc[:,:12].sort_values('mu_KD', ascending=False)

# SCR
scr_specific = scr.merge(kd, how='left', on=['chrI', 'startI', 'endI', 'chrJ', 'startJ', 'endJ'], suffixes=['_SCR', '_KD'])
test = scr_specific.loc[:,scr_specific.columns.str.contains('KD')].isna().sum(axis=1) != 0
scr_specific = scr_specific.loc[test,:].iloc[:,:12].sort_values('mu_SCR', ascending=False)

# Common
common = kd.merge(scr, how='inner', on=['chrI', 'startI', 'endI', 'chrJ', 'startJ', 'endJ'], suffixes=['_KD', '_SCR'])
kd_specific.shape[0] + common.shape[0] + scr_specific.shape[0]


## 


# SCR specific and CtIP
enh = pd.read_csv(path_enhancers, sep='\t', header=None)
enh.columns = ['chrI', 'startI', 'endI', 'cluster']
enh = enh.dropna()

# Filter out CtIP enancers
df_coords = filter_enhancers(enh, scr_specific)

# Intersect and rank mutated enhancers
df_muts = pd.read_csv(os.path.join(path_results, 'andrea.CtIP_mutated_enhancers.tsv'), sep='\t')


a = df_coords['chr'] + ':' + df_coords['enhancer_summit'].astype(np.int64).astype('str')
b = df_muts['seqnames'] + ':' + df_muts['summit'].astype('str')


a.isin(b)

a


