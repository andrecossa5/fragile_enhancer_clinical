"""
Utils for survival analysis.
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxPHSurvivalAnalysis


##


def _univariate_Cox(df, gene=None, status_cov='RFS_STATUS', time_cov='RFS_time'):
    """
    Perform simple univariate analysis.
    """

    # Get X and y
    test = ~df[gene].isna()
    X = df.loc[test, gene].values[:, np.newaxis]                                # Remove NAs, get gene column, to 2D
    X = StandardScaler().fit_transform(X)                                       # Scale gene column
    y = np.zeros((X.shape[0],), dtype=[('delta', '?'), ('time', 'f8')]) 
    y[['delta']] = np.where(df.loc[test, status_cov], 1, 0)
    t = df.loc[test, time_cov].values
    y[['time']] = (t - t.mean()) / t.std()                                      # Scale time column

    # Here we go
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=1234)

    hratio = []
    c_index = []
    for train_idx, test_idx in splitter.split(np.arange(y.size), y['delta']):
        model = CoxPHSurvivalAnalysis(n_iter=1000)
        model.fit(X[train_idx,:], y[train_idx])
        h = np.exp(model.coef_)
        c = model.score(X[test_idx,:], y[test_idx])
        hratio.append(h)
        c_index.append(c)

    return np.mean(hratio), np.std(hratio), np.mean(c_index)


##


def _multivariate_Cox(df, gene=None, status_cov='RFS_STATUS', time_cov='RFS_time', confounders=None):
    """
    Perform multivariate analysis with confounders.
    """

    # Get X and y
    cols = [gene] + confounders
    idx = df.loc[:, cols].dropna().index                        
    X = StandardScaler().fit_transform(df.loc[idx, cols].values)                # Scale gene column
    y = np.zeros((X.shape[0],), dtype=[('delta', '?'), ('time', 'f8')]) 
    y[['delta']] = np.where(df.loc[idx, status_cov], 1, 0)
    t = df.loc[idx, time_cov].values
    y[['time']] = (t - t.mean()) / t.std()                                      # Scale time column

    # Here we go
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=1234)

    hratio = []
    c_index = []
    for train_idx, test_idx in splitter.split(np.arange(y.size), y['delta']):
        model = CoxPHSurvivalAnalysis(n_iter=1000)
        model.fit(X[train_idx,:], y[train_idx])
        h = np.exp(model.coef_)[:,np.newaxis]
        c = model.score(X[test_idx,:], y[test_idx])
        hratio.append(h)
        c_index.append(c)

    l_ = [ f'CV{i}' for i in range(1,6) ]
    hratios = pd.DataFrame(np.concatenate(hratio, axis=1), index=cols, columns=l_)
    hratios = pd.concat([hratios.mean(axis=1), hratios.std(axis=1)], axis=1)
    hratios.columns = ['mean', 'std']
    c_index = np.mean(c_index)

    return hratios, c_index


##