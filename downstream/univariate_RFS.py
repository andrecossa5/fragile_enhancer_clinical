#Univariate analysis OS, within METABRIC data

########################################################################

##Standard modules and utilities
import copy
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from sksurv.preprocessing import OneHotEncoder
from itertools import product

##Models
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest

##Model selection and evaluation
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.compare import compare_survival

#My utilities
import sys
sys.path.append('/Users/IEO5505/Desktop/ML_TDA/S_3/scripts/')
import utilities as my

#NO WARNINGs
import warnings
warnings.filterwarnings("ignore")

########################################################################

#All METABRIC

##Load data
df = pd.read_csv('/Users/IEO5505/Desktop/ML_TDA/S_3/data/new/final_df_subsetted.csv', sep=',')

##Prep data
X, y, X_train, X_test, y_train, y_test, features = my.prepSingle_sv(df, 'RFS', 'METABRIC', only_genes = False)
clinical = features[:8]
genes = features[8:]

##Univariate analysis 
folds = 5
C_score, hratio = my.univariate_Cox(X, y, features, folds)

##Save
C_score.sort_values(by = 'mean', ascending = False).to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/C_score_univariate_RFS_all.xlsx')
hratio.sort_values(by = 'mean', ascending = False).to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/hratio_univariate_RFS_all.xlsx')
hratio = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/hratio_univariate_RFS_all.xlsx', index_col=0)
C_score = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/C_score_univariate_RFS_all.xlsx', index_col=0)

##Inspect results
###C_score: all, clinical, genes
C_score.sort_values(by = 'mean', ascending = False).head(n = 20)
C_score.loc[clinical, :].sort_values(by = 'mean', ascending = False).head(n = 20)
C_score.loc[genes, :].sort_values(by = 'mean', ascending = False).head(n = 20)

###hratio: all, clinical, genes
hratio.sort_values(by = 'mean', ascending = False).head(n = 20)
hratio.loc[clinical, :].sort_values(by = 'mean', ascending = False).head(n = 20)
hratio.loc[genes, :].sort_values(by = 'mean', ascending = False).head(n = 20)

##hratio and C-score of all features and genes
sum(C_score['mean'] > 0.5)
sum(C_score.loc[genes, :]['mean'] > 0.5)
sum(hratio['mean'] > 1.0)
sum(hratio.loc[genes, :]['mean'] > 1.0)

##Save 'good' genes
by_c = list(C_score.loc[genes][C_score['mean'] > 0.5].index)
by_hratio = list(hratio.loc[genes][hratio['mean'] > 1.0].index)
all_genes = pd.DataFrame({'genes' : list(set(by_c) & set(by_hratio))})
all_genes['mean_C_score'] = list(C_score.loc[list(all_genes['genes'])]['mean'])
all_genes['mean_C_hratio'] = list(hratio.loc[list(all_genes['genes']), 'mean'])
all_genes.to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/candidates_genes_all.xlsx')


##


##############################------------------------------------------

#Only basal (only genes)

##Load data
df = pd.read_csv('/Users/IEO5505/Desktop/ML_TDA/S_3/data/new/final_df_subsetted_basal.csv', sep=',')

##Prep data (only genes here)
X, y, X_train, X_test, y_train, y_test, features = my.prepSingle_sv(df, 'RFS', 'METABRIC', only_genes = True)
features

##Univariate analysis 
folds = 5
C_score, hratio = my.univariate_Cox(X, y, features, folds)

##Inspect results
C_score.sort_values(by = 'mean', ascending = False).head(n = 20)
hratio.sort_values(by = 'mean', ascending = False).head(n = 20)

##Save
C_score.sort_values(by = 'mean', ascending = False).to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/C_score_univariate_RFS_basal.xlsx')
hratio.sort_values(by = 'mean', ascending = False).to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/hratio_univariate_RFS_basal.xlsx')
hratio = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/hratio_univariate_RFS_basal.xlsx', index_col=0)
C_score = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/C_score_univariate_RFS_basal.xlsx', index_col=0)

##hratio and C-score of all features
sum(hratio['mean'] > 1.0)
sum(C_score['mean'] > 0.5)

##Save 'good' genes
by_c = list(C_score.loc[features][C_score['mean'] > 0.5].index)
by_hratio = list(hratio.loc[features][hratio['mean'] > 1.0].index)
basal_genes = pd.DataFrame({'genes' : list(set(by_c) & set(by_hratio))})
basal_genes['mean_C_score'] = list(C_score.loc[list(basal_genes['genes'])]['mean'])
basal_genes['mean_C_hratio'] = list(hratio.loc[list(basal_genes['genes']), 'mean'])
basal_genes.to_excel('/Users/IEO5505/Desktop/ML_TDA/S_3/plots_results/univariate_RFS/candidates_genes_basal.xlsx')


##


'''
NBBB: Summary for this analysis.

All 
-132 feats (128 genes over 221) have hr > 1; 164 feats (157 genes over 231) have c-index >= 0.5

Basal
-138 over 221 genes have hr > 1; 117 over 231 have c-index >= 0.5
'''


########################################################################




