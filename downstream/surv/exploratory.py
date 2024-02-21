"""
Exploratory data analysis. How much do we have??
"""

import os
import numpy as np
import pandas as pd


##


# Paths
path_main = '/Users/IEO5505/Desktop/fragile-enchancer-clinical'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'exploratory')

# Subpaths
path_clinical = os.path.join(path_data, 'clinical')
path_WES = os.path.join(path_data, 'WES')
path_WGS = os.path.join(path_data, 'WGS')


##


#1 How many patients, specimens and samples? How many with clinical annot?
donor_df = pd.read_csv(os.path.join(path_clinical, 'donor.csv'), index_col=0)
sample_df = pd.read_csv(os.path.join(path_clinical, 'sample.csv'), index_col=0)
specimen_df = pd.read_csv(os.path.join(path_clinical, 'specimen.csv'), index_col=0)


##


# All numbers
donor_df.index.unique().size        # 789 
sample_df.index.unique().size       # 3339
specimen_df.index.unique().size     # 1592

# Specimen types
specimen_types_all = specimen_df['specimen_type'].value_counts()
specimen_types_all.sum()
specimen_types_all.to_csv(os.path.join(path_results, 'specimen_type.csv'))

test = (~donor_df['donor_vital_status'].isna()) & \
       (~donor_df['donor_age_at_last_followup'].isna()) & \
       (~donor_df['donor_age_at_diagnosis'].isna())
donors_with_clinical = donor_df.loc[test].index
donors_with_clinical.size                             # 176

specimen_type_with_clinical = (
    specimen_df
    .query('icgc_donor_id in @donors_with_clinical')
    ['specimen_type'].value_counts()
)
specimen_type_with_clinical.sum()                    # 366
specimen_type_with_clinical.to_csv(
    os.path.join(path_results, 'specimen_type_with_clinical_info.csv')
)


##


# 2 How many patients, specimens and samples with SNPs, SVs and indels?
# 3 How many patients, specimens and samples with SNPs, SVs and indels and clinical info?
SNPs = pd.read_csv(os.path.join(path_WGS, 'SNPs.csv'), index_col=0)
INDELs = pd.read_csv(os.path.join(path_WGS, 'INDELs.csv'), index_col=0)
SVs = pd.read_csv(os.path.join(path_WGS, 'SVs.csv'), index_col=0)

pd.Series(SNPs['icgc_donor_id'].unique()).isin(donor_df.index).sum()            # 191
pd.Series(SNPs['icgc_donor_id'].unique()).isin(donors_with_clinical).sum()      # 68
(
    SNPs.query('icgc_donor_id in @donors_with_clinical')
    .groupby('cluster')
    .size().sort_values(ascending=False)
    .to_csv(os.path.join(
        path_results, 'n_SNPs_per_enh_cluster_in_donors_with_clinical.csv'
    ))
)

pd.Series(INDELs['icgc_donor_id'].unique()).isin(donor_df.index).sum()          # 123
pd.Series(INDELs['icgc_donor_id'].unique()).isin(donors_with_clinical).sum()    # 10 
(
    INDELs.query('icgc_donor_id in @donors_with_clinical')
    .groupby('cluster')
    .size().sort_values(ascending=False)
    .to_csv(os.path.join(
        path_results, 'n_INDELs_per_enh_cluster_in_donors_with_clinical.csv'
    ))
)

pd.Series(SVs['icgc_donor_id'].unique()).isin(donor_df.index).sum()             # 412
pd.Series(SVs['icgc_donor_id'].unique()).isin(donors_with_clinical).sum()       # 0


##



