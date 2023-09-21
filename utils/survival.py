"""
Utils for survival analysis.
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxPHSurvivalAnalysis


##


def one_hot_from_labels(df, var_name):
    """
    My one_hot encoder from a categorical-like pd.Series.
    """
    y = df[var_name]
    labels = y.unique()
    labels = [ x for x in labels if x != 'unknown' ]

    if len(labels)>2:
        Y = np.concatenate(
            [ np.where(y == x, 1, 0)[:, np.newaxis] for x in labels ],
            axis=1
        )
        one_hot_df = pd.DataFrame(
            Y, index=df.index, columns=[f'{var_name}_{x}' for x in labels]
        )
    else:
        Y = np.where(y == labels[0], 1, 0)
        one_hot_df = pd.Series(Y).to_frame(var_name)
    
    return one_hot_df


##


def univariate_Cox(df, molecular_feature=None, status_cov='RFS_STATUS', time_cov='RFS_time'):
    """
    Perform simple univariate analysis.
    """

    # Get X and y
    test = ~df[molecular_feature].isna()
    X = df.loc[test, molecular_feature].values[:, np.newaxis]               # Remove NAs, get gene column, to 2D
    X = scale(X)                                                            # Scale gene column
    y = np.zeros((X.shape[0],), dtype=[('delta', '?'), ('time', 'f8')]) 
    y[['delta']] = np.where(df.loc[test, status_cov], 1, 0)
    y[['time']] = df.loc[test, time_cov].values                                

    # Here we go
    splitter = StratifiedKFold(n_splits=5, shuffle=True)                    # No seed??

    hratio = []
    c_index = []
    for train_idx, test_idx in splitter.split(np.arange(y.size), y['delta']):
        model = CoxPHSurvivalAnalysis(n_iter=1000)
        model.fit(X[train_idx,:], y[train_idx])
        h = np.exp(model.coef_[0])
        c = model.score(X[test_idx,:], y[test_idx])
        hratio.append(h)
        c_index.append(c)

    s = pd.Series({
        'model_type' : 'univariate',
        'molecular_feature' : molecular_feature,
        'HR_mean' : np.mean(hratio),
        'HR_std' : np.std(hratio),
        'C_mean' : np.mean(c_index),
        'C_std' : np.std(c_index),
        'HR_confounders' : None
    })

    return s


##


def multivariate_Cox(df, molecular_feature=None, status_cov='RFS_STATUS', 
                    time_cov='RFS_time', confounders=None):
    """
    Perform multivariate analysis with confounders.
    """

    # X
    df_ = df[[molecular_feature]+confounders].dropna() 
    for x in confounders:
        s = df_[x]
        if pd.api.types.is_numeric_dtype(s):
            df_.loc[:,x] = scale(s)
        else:
            df_ = df_.drop(x, axis=1).join(one_hot_from_labels(df_, x)) 
    var_names = df_.columns   
    X = df_.values

    # y
    df_ = df.loc[df_.index]
    y = np.zeros((df_.shape[0],), dtype=[('delta', '?'), ('time', 'f8')]) 
    y[['delta']] = np.where(df_[status_cov], 1, 0)
    y[['time']] = df_[time_cov].values                               # Scale time column?? 
    
    # Here we go
    splitter = StratifiedKFold(n_splits=5, shuffle=True)

    hratio = []
    c_index = []
    for train_idx, test_idx in splitter.split(np.arange(y.size), y['delta']):
        model = CoxPHSurvivalAnalysis(n_iter=1000)
        model.fit(X[train_idx,:], y[train_idx])
        h = np.exp(model.coef_)
        c = model.score(X[test_idx,:], y[test_idx])
        hratio.append(h)
        c_index.append(c)

    df_confounders = pd.DataFrame(
        np.concatenate([ x[1:][:,np.newaxis] for x in hratio], axis=1),
        index=var_names[1:]
    )
    df_confounders = df_confounders.T.describe().loc[['mean', 'std']]
    df_confounders.index = df_confounders.index.map(lambda x: f'HR_{x}')

    s = pd.Series({
        'model_type' : 'multivariate',
        'molecular_feature' : molecular_feature,
        'HR_mean' : np.mean([ x[0] for x in hratio]),
        'HR_std' : np.std([ x[0] for x in hratio]),
        'C_mean' : np.mean(c_index),
        'C_std' : np.std(c_index),
        'HR_confounders' : df_confounders.to_dict()
    })

    return s


##