"""
Script for Metabric data wrangling.
"""

import os
import numpy as np
import pandas as pd


##


# Paths
path_data = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/METABRIC/raw'
path_results = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/METABRIC'


##


# Read and basic cleaning
patients_df = pd.read_csv(
    os.path.join(path_data, 'data_clinical_patient.txt'),   # patients
    index_col=0, sep='\t'
) 
to_retain = ["AGE_AT_DIAGNOSIS", "OS_STATUS", "OS_MONTHS", "RFS_STATUS", "RFS_MONTHS"]
patients_df = patients_df.loc[:, to_retain].dropna().copy()

samples_df = pd.read_csv(
    os.path.join(path_data, 'data_clinical_sample.txt'),    # samples
    index_col=0, sep='\t'
)
to_retain = ['ER_STATUS', 'HER2_STATUS', 'PR_STATUS', 'TUMOR_STAGE']
samples_df = samples_df.loc[:, to_retain].dropna().copy()

M = pd.read_csv(
    os.path.join(path_data, 'data_expression_median.txt'),   # expression
    index_col=0, sep='\t'
)
to_drop = M.index.value_counts().loc[lambda x: x>1].index.tolist()
M = M.drop(index=to_drop, columns=['Entrez_Gene_Id']).T.copy()
M = M.astype(np.float16)


##


# Join
df = patients_df.join([samples_df, M])

# Filter and format colnames
df.loc[df['OS_STATUS'] == '0:LIVING', 'OS_STATUS'] = False
df.loc[df['OS_STATUS'] == '1:DECEASED', 'OS_STATUS'] = True
df.loc[df['RFS_STATUS'] == '0:Not Recurred', 'RFS_STATUS'] = False
df.loc[df['RFS_STATUS'] == '1:Recurred', 'RFS_STATUS'] = True
df.insert(df.columns.get_loc('RFS_STATUS')+1, 'OS_time', df['OS_MONTHS']*31)
df.insert(df.columns.get_loc('OS_time')+1, 'RFS_time', df['RFS_MONTHS']*31)
df = df.loc[:,~df.columns.str.contains('MONTHS')].copy()

_L = [
    (df['ER_STATUS'] == 'Negative') & (df['HER2_STATUS'] == 'Negative') \
    & (df['PR_STATUS'] == 'Negative'),
    (df['ER_STATUS'] == 'Negative') & (df['HER2_STATUS'] == 'Positive') \
    & (df['PR_STATUS'] == 'Negative'),
    (df['ER_STATUS'] == 'Positive') & (df['HER2_STATUS'] == 'Negative') \
    & (df['PR_STATUS'] == 'Positive'),
]
s = np.select(_L, ['TNBC', 'Her2', 'Luminal'], default='unknown')
df.insert(df.columns.get_loc('TUMOR_STAGE')+1, 'subtype', s)


##


# Save
df.to_csv(os.path.join(path_results, 'formatted_data.csv'))

