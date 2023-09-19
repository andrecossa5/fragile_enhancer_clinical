"""
Script for univariate survival analysis with DE results.
"""

import os
from .utils.survival import *


##


# Paths
path_data = '/Users/IEO5505/Desktop/fragile-enchancer-clinical/data/'

# Read formatted METABRIC data and DEG list 
df = pd.read_csv(os.path.join(path_data, 'METABRIC', 'formatted_data.csv'), index_col=0)
DEGs = pd.read_csv(os.path.join(path_data, 'expression', 'GRHL2_KD_DEG.csv'), index_col=0)
DEGs = DEGs.loc[DEGs['FDR']<.1,['logFC', 'FDR']].sort_values('FDR').copy() # Remove FDR >= .1
clinical_features = df.columns[:10].tolist()
expression_features =  df.columns[df.columns.isin(DEGs.index)].tolist()


##


# Perform univariate analysis
gene = expression_features[-1]
status_cov = 'RFS_STATUS'
time_cov  = 'RFS_time'
confounders = ['AGE_AT_DIAGNOSIS', 'TUMOR_STAGE']


##


# Calculate survival indeces

# Univariate
d = {
    x : _univariate_Cox(df, gene=x, status_cov='RFS_STATUS', time_cov='RFS_time') \
    for x in expression_features[:10]
}
df_results = pd.DataFrame(d, index=['HR_mean', 'HR_std', 'C_mean']).T

# Multivariate
d = {
    x : _multivariate_Cox(df, gene=x, status_cov='RFS_STATUS', time_cov='RFS_time', confounders=['AGE_AT_DIAGNOSIS']) \
    for x in expression_features[:10]
}

# Save results



##


# Load and visualize with KMs.






















## 