"""
Get the set of loops from SCR and KD.
"""

import os
import pandas as pd
import numpy as np


## 


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical'
path_data = os.path.join(path_main, 'data', 'Hi_Chip', 'loops')
path_results = os.path.join(path_main, 'results', 'Hi_Chip')

# Read interactions
SCR = pd.read_csv(os.path.join(path_data, 'SCR_loops.bed'), sep='\t')

SCR.query('chr1==11')['e2'].describe()



# SCR = pd.read_csv(os.path.join(path_data, 'SCR_loops_8k.bed'), sep='\t')
KD = pd.read_csv(os.path.join(path_data, 'KD_loops.bed'), sep='\t')
enh = pd.read_csv(os.path.join(path_main, 'data', 'Hi_Chip', 'CtIP_enhancers.txt'), sep='\t', header=None)
enh.columns = ['chr', 'start', 'stop', 'cluster']
enh = enh.query('chr!="chrX" and chr!="#"')             # NB
enh['chr'] = enh['chr'].str.strip('chr').astype(int)

# Reformat
SCR['loop'] = 'chr' + SCR['chr1'].astype('str') + '_' + SCR['s1'].astype('str') + '_' \
               + SCR['e1'].astype('str') + '_chr_' + SCR['chr2'].astype('str') + '_'  \
               + SCR['s2'].astype('str') + '_' + SCR['e2'].astype('str')
KD['loop'] = 'chr' + KD['chr1'].astype('str') + '_' + KD['s1'].astype('str') + '_' \
               + KD['e1'].astype('str') + '_chr_' + KD['chr2'].astype('str') + '_'  \
               + KD['s2'].astype('str') + '_' + KD['e2'].astype('str')

# Relationship
sign_test = lambda x: x['Q-Value_Bias']<=0.1
SCR.loc[sign_test].shape[0]
KD.loc[sign_test].shape[0]
loops_SCR = set(SCR.loc[sign_test]['loop'].values)
loops_KD = set(KD.loc[sign_test]['loop'].values)
len(loops_SCR & loops_KD)
len(loops_SCR - loops_KD)
len(loops_KD - loops_SCR)


##


# Enhancer contacts
# enh_loops_L = []
# for i in range(enh.shape[0]):
#     s = enh.iloc[i,:]
#     start = s['start']
#     chrom = s['chr']
#     kd = SCR.loc[sign_test].query('chr1==@chrom and chr2==@chrom')
#     cont1 = (start>=kd['s1']) & (start<=kd['e1'])
#     cont2 = (start>=kd['s2']) & (start<=kd['e2'])
#     n = np.sum(cont1 | cont2)
#     if n>0:
#         if np.sum(cont1)>0:
#             df_ = kd.iloc[np.where(cont1)[0],:]
#         elif np.sum(cont2)>0:
#             df_ = kd.iloc[np.where(cont2)[0],:]
#         enh_loops_L.append(df_.assign(enhancer=f'enhancer_chr{chrom}_{round(start)}', enh_cluster=s['cluster']))
# df_loops = pd.concat(enh_loops_L)

# Save
# df = df_loops[['loop', 'enhancer', 'enh_cluster', 'Coverage1', 'Coverage2', 'cc', 'Bias1', 'Bias2', 'Q-Value_Bias']]
# df.to_csv(os.path.join(path_results, 'SCR_CtIP_enh_loops_8k.csv'))


##


# Read filtered loops
KD_loops = pd.read_csv(os.path.join(path_results, 'KD_CtIP_enh_loops.csv'), index_col=0)
SCR_loops = pd.read_csv(os.path.join(path_results, 'SCR_CtIP_enh_loops_8k.csv'), index_col=0)

KD_loops.shape
SCR_loops.shape

common = list(set(KD_loops['loop'].to_list()) & set(SCR_loops['loop'].to_list()))
new = list(set(KD_loops['loop'].to_list()) - set(SCR_loops['loop'].to_list()))
depleted = list(set(SCR_loops['loop'].to_list()) - set(KD_loops['loop'].to_list()))


SCR_loops.set_index('loop').loc[common]
SCR


##
