"""
Integrated analysis of HiCHip loops.
"""

import os
from utils.integrated import *
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical'
path_func_genomics = os.path.join(path_main, 'data', 'functional_genomics')
path_hichip = os.path.join(path_func_genomics, 'HiChip')
path_results = os.path.join(path_main, 'results', 'integrated')
path_enhancers = os.path.join(path_main, 'data', 'functional_genomics', 'others', 'CtIP_enhancers.txt')
path_enhancers = os.path.join(path_main, 'data', 'functional_genomics', 'others', 'Cluster_GRHL_Enh_All.txt')
path_tss = os.path.join(path_func_genomics, 'others', 'TSSs_elisa', 'TSSs_from_hg19_ncbiRefSeq.with_activity.tsv')
path_expression = os.path.join(path_main, 'data', 'expression', 'KD_vs_SCR_DEGs.csv')
res = 4


##


# Read loops
scr = pd.read_csv(os.path.join(path_hichip, 'filtered_loops', 'SCR', f'SCR_loops_{res}_kb.tsv.gz'), sep='\t')
kd = pd.read_csv(os.path.join(path_hichip, 'filtered_loops', 'KD', f'KD_loops_{res}_kb.tsv.gz'), sep='\t')

# KD
kd_specific = kd.merge(
    scr, on=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'], 
    how='left', suffixes=['_KD', '_SCR']
)
test = kd_specific.loc[:,kd_specific.columns.str.contains('SCR')].isna().sum(axis=1) != 0
kd_specific = kd_specific.loc[test,:].iloc[:,:12].sort_values('mu_KD', ascending=False)

# SCR
scr_specific = scr.merge(
    kd, on=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'], 
    how='left', suffixes=['_SCR', '_KD']
)
test = scr_specific.loc[:,scr_specific.columns.str.contains('KD')].isna().sum(axis=1) != 0
scr_specific = scr_specific.loc[test,:].iloc[:,:12].sort_values('mu_SCR', ascending=False)

# Common
common = kd.merge(
    scr, on=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'], 
    how='inner', suffixes=['_KD', '_SCR']
)
kd_specific.shape[0] + common.shape[0] + scr_specific.shape[0]


## 


# SCR specific loops: CtIP enhancers and other interactors
enh = pd.read_csv(path_enhancers, sep='\t', header=None)
enh.columns = ['chr', 'start', 'end', 'cluster']
enh = enh.dropna()
enh['cluster'].value_counts()

# SCR specific (KO depleted) loops which connects CtIP enhancers with other interactors (ctip_loops)
ctip_loops = filter_loops_with_enhancers(scr_specific, enh)
ctip_loops.shape[0] / scr_specific.shape[0]
# ctip_loops.to_csv(os.path.join(path_results, 'all_GRHL2_loops.tsv'), sep='\t', index=False)

# ctip_loops CtIP cluster prevalence
ctip_loops['isin_bin'].value_counts()
ctip_loops['cluster'].value_counts(normalize=True)

# ctip_loops connecting ctip-enh with known gene tss 
df_tss = pd.read_csv(path_tss, sep='\t').query('activity=="active"')
ctip_tss_loops = (
    scr_specific.merge(
        filter_loops_with_tss(ctip_loops, df_tss),
        how='inner', on=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'] 
    )
    .sort_values('mu_SCR', ascending=False)
)
# ctip_tss_loops.to_csv(os.path.join(path_results, 'GRHL2_active_tss_loops.tsv'), sep='\t', index=False)

# How many ctip-enh --> tss loops, per gene?
ctip_tss_loops['gene'].value_counts().describe()
ctip_tss_loops['gene'].value_counts().head(50)
ctip_tss_loops['gene'].values

# DEGs
df_de = pd.read_csv(path_expression)
ctip_tss_de_loops = ctip_tss_loops.loc[ctip_tss_loops['gene'].isin(df_de['GeneName'])]
# ctip_tss_de_loops.to_csv(os.path.join(path_results, 'GRHL2_active_tss_DE_loops.tsv'), sep='\t', index=False)
ctip_tss_de_loops.shape[0]
g1 = set(ctip_tss_loops['gene'].to_list())
g2 = set(df_de['GeneName'].to_list())
genes = list(g1 & g2)
df_de.loc[df_de['GeneName'].isin(genes)][['GeneName','logFC', 'logCPM']].sort_values('logFC')['GeneName'].to_list()

# Save some loops
ctip_tss_de_loops.to_csv(os.path.join(path_results, 'prova_loops_4kb.csv'))
g1 = set(pd.read_csv(os.path.join(path_results, 'prova_loops_4kb.csv'), index_col=0)['gene'].unique().tolist())
g2 = set(pd.read_csv(os.path.join(path_results, 'prova_loops_8kb.csv'), index_col=0)['gene'].unique().tolist())
len(g1&g2)


##


# Intersect and rank mutated enhancers
# df_muts = pd.read_csv(os.path.join(path_results, 'andrea.CtIP_mutated_enhancers.tsv'), sep='\t')
# 
# 
# a = df_coords['chr'] + ':' + df_coords['enhancer_summit'].astype(np.int64).astype('str')
# b = df_muts['seqnames'] + ':' + df_muts['summit'].astype('str')
# 
# 
# a.isin(b)


