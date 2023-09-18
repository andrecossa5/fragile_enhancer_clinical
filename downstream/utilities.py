#Useful functions for ML_pipelines

#####################################################################

#Import modules and functions

##Standard modules and utilities
import copy
import random
import multiprocess as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from sksurv.preprocessing import OneHotEncoder
from itertools import product

##Models
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, ComponentwiseGradientBoostingSurvivalAnalysis, GradientBoostingSurvivalAnalysis
from lifelines import CoxPHFitter

##Model selection and evaluation
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.compare import compare_survival

######################################################################

#Define utilities:

#Survival analysis functions

########Example data
df = pd.read_csv('/Users/IEO5505/Desktop/ML_TDA/Survival_Analysis/TCGA_and_METABRIC/data/new/final_df_subsetted.csv', sep=',')
########


##


##prepSingle_sv(): prepare survival analysis ingredients from a single database;
def prepSingle_sv(df, kind, database, only_genes = False):
    '''
    Prep X and y data manipulating the original data table, and perform train-test splitting.
    '''

    ##Get rid of identifiers, and surivial columns not specific of the analysis kind
    if kind == 'OS':
        df = df.drop(['sample', 'case', 'Prog_free_status', 'Prog_free_time'], axis = 1)
        df = df.rename(columns = {'OS_status' : 'status', 'OS_time' : 'time'})
    elif kind == 'RFS':
        df = df.drop(['sample', 'case', 'OS_status', 'OS_time'], axis = 1)
        df = df.rename(columns = {'Prog_free_status' : 'status', 'Prog_free_time' : 'time'})
    
    ##Take out X and y of features and targets. NB: X is a standard np.array whose clinical 
    ##features have to be encoded in OneHot scheme, while y is a structured array, 
    ##storing an indicator and a time, depending on the analysis 'kind'.

    ###X 
    X = df.drop(['status', 'time'], axis = 1)
    
    ###Feature choice
    if database == 'TCGA':
        if only_genes == False: 
            for col in ['stage', 'subtype', 'ajcc_m_score', 'race']:
                X[col] = df[col].astype('category')
            X = OneHotEncoder().fit_transform(X)

        else:
            X = X.drop(['age_at_diagnosis', 'stage', 'subtype', 'ajcc_m_score', 'race'], axis=1) 
    elif database == 'METABRIC':
        if only_genes == False: 
            for col in ['stage', 'subtype']: 
                X[col] = df[col].astype('category')
            X = OneHotEncoder().fit_transform(X)
    
        else:
            X = X.drop(['age_at_diagnosis', 'stage', 'subtype'], axis=1) 

    ###Store info about encoded features
    features = list(X.columns)

    ###Finally convert X to numpy array
    X = X.values

    ###y
    y = np.zeros((X.shape[0],), dtype = [('delta', '?'), ('time', 'f8')]) 
    y[['delta']] = df.status.values
    y[['time']] = df.time.values

    ##Train-test split. StratifiedShuffleSplit: 2 splits, 80-20% between training and validation, 
    ##each time
    sss = StratifiedShuffleSplit(n_splits = 2, test_size = 0.2, random_state = 0)
    ##Index X and y, test and train
    for train_index, test_index in sss.split(X = np.zeros(len(y['delta'])), y = y['delta']):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

    return(X, y, X_train, X_test, y_train, y_test, features)
 

##


######Test
X, y, X_train, X_test, y_train, y_test, features = prepSingle_sv(df, 'OS', 'METABRIC', only_genes = False)
######


##


########Example
folds = 5
X.shape
len(features)
#########

#univariate_Cox(): given the list of features, this function fit a univariate Cox models for each
#each feature and gives CV results 
def univariate_Cox(X_train, y_train, features, folds = 5):
    '''
    Given X and y, the function fit an univariate CoxPH model for each feature of the training set, 
    using stratified 5 fold CV. Each time, the model is evaluated using C-index. The mean and std over CV
    and the mean log-hazard ratio are reported.
    ''' 

    ##Initialize C_score and hratio arrays
    C_score = np.full((len(features), folds), 1.111)
    hratio = np.full((len(features), folds), 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over features and folds indeces. Fit on tr; score model on t; and store infos.
    for i in range(len(features)): 
        fold_c = 0
        for tr, t in skf.split(np.zeros(len(y_train)), y_train['delta']):
            m = CoxPHSurvivalAnalysis(n_iter = 1000)
            m.fit(X_train[tr, i:i+1], y_train[tr])
            hratio[i, fold_c] =  np.exp(m.coef_)
            C_score[i, fold_c] = m.score(X_train[t, i:i+1], y_train[t])
            fold_c += 1
        print('Done for feature %s' % features[i])
    
    #Convert to df and calculate mean and std
    C_score = pd.DataFrame(C_score, index = features, columns = [ 'fold_' + str(i) for i in range(5) ])
    C_score['mean'] = C_score.mean(axis=1)
    C_score['std'] = C_score.std(axis=1)
    hratio = pd.DataFrame(hratio, index = features, columns = [ 'fold_' + str(i) for i in range(5) ])
    hratio['mean'] = hratio.mean(axis=1)
    hratio['std'] = hratio.std(axis=1)

    return (C_score, hratio)


########Test
#C_score, hratio = univariate_Cox(X, y, features, folds = 5)
#######


##


########Example
#Good features
all_genes = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/TCGA_and_METABRIC/plots_results/univariate_OS/candidates_genes_all.xlsx', index_col=0)
#basal = pd.read_excel('/Users/IEO5505/Desktop/ML_TDA/S_2/plots_results/TCGA/univariate_OS_TCGA/candidates_genes_basal.xlsx', index_col=0)
####basal['p'] = (2 * basal.mean_C_score * basal.mean_C_hratio) / (basal.mean_C_score + basal.mean_C_hratio)
####set(basal.sort_values('mean_C_score', ascending=False)['genes'][:7]) & set(basal.sort_values('mean_C_hratio', ascending=False)['genes'][:7])
####basal[basal['genes'] == 'YPEL3']
#all_genes.head()
all_genes = list(all_genes['genes'])
#basal = list(basal['genes'])
##seed_all = 'C8orf33'
####seed_basal = 'YPEL3'
####Set of indeces for X subsetting
#idx_features_to_use_basal = [features.index(x) for x in basal]
idx_features_to_use_all = [features.index(x) for x in all_genes]
len(idx_features_to_use_all)
####Subset X
#X_basal = X_train[:, idx_features_to_use_basal]
X_all = X_train[:, idx_features_to_use_all]
###Seed index in X_all
#seed = [features[i] for i in idx_features_to_use_all].index(seed_all)
########


##


#internal_evaluator(): implements evalutation of a candidate 'signature'.
def internal_evaluator(a, y, folds):
    '''
    Given an array a storing k covariates values, the function evaluate by 5-fold CV the 
    predictive power of such a signature, computing the C_index of the resulting CoxPH
    models.
    '''  
    
    ##Initialize C_score and hratio lists
    C = np.full(folds, 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over folds indeces. Fit on tr; score model on t; and store infos.
    fold_c = 0
    for tr, t in skf.split(np.zeros(len(y)), y['delta']):
        m = CoxPHSurvivalAnalysis(n_iter = 1000000)
        m.fit(a[tr, ], y[tr])
        C[fold_c] = m.score(a[t, ], y[t])
        fold_c += 1
    C_mean = np.mean(C)

    return(C_mean)

   
##


#make_combos(): make all the combos at k + 1, using the remaining indeces.
def make_combos(signature, n):
    '''
    Given an existing signature and some other features, produces a list of new combinations
    to test.
    ''' 

    #Take all the features indeces that are not in signature, and append it to it 
    #one at the time, to produce a list of new combos (lists) to test.
    others = [[i] for i in list(range(n)) if i not in signature]
    new_combos = []
    for x in others:
        new_combos.append(signature + x)

    return(new_combos)

    
##


#forward_selector(): implements the first strategy for signature builiding: forward selection.
def forward_selector(X, y, seed, folds, path):
    '''
    Given X, y training data the function first fit a univariate model for the seed feature.
    Then, it adds iteratively each of the remaining features to the predictor, and test each model 
    C_index with 5 fold CV. The best model C_score is then tested >= to previous k-1 step best one,
    and if the test results positive, another iteration is allowed, otherwise the 'optimal' list of
    predictors is returned. 
    '''

    #Initialize: 1) results: combo_tested, k, C_index dictionary 2) combo_chosen.
    results = {}
    combo_chosen = []

    #Test initial signature (k = 1, seed gene), store its C_index
    k = 1
    results['combo_tested'] = [[seed]]
    results['k'] = [k]
    a = X[:, seed:seed+1]
    results['C_index'] = [internal_evaluator(a, y, folds)]
    #Append the seed combo to combo_chosen
    combo_chosen.append([seed])
    
    #Iteratively, until the best new has better C than the previous one 
    proceed = True
    n = X.shape[1]
    while proceed:
        #1) Create a list of new combo (of k + 1 genes)
        new_combos = make_combos(combo_chosen[k-1], n)
        #2) For each new combo --> CV testing, scoring. Append outputs to results
        k += 1
        for combo in new_combos:
            results['combo_tested'].append(combo)
            results['k'].append(k)
            a = X[:, combo]
            results['C_index'].append(internal_evaluator(a, y, folds))
        #3) Choose the best combo among the last k combos
        look_at = [ i for i in list(range(len(results['k']))) if results['k'][i] == k ]
        new_combos_results = [ results['C_index'][i] for i in look_at ]
        best_new_c = max(np.array(new_combos_results))
        best_new = new_combos[ new_combos_results.index(best_new_c) ]
        #4) Decide to proceed or not
        proceed = best_new_c > results['C_index'][ results['combo_tested'].index(combo_chosen[k-2]) ]
        #if yes, add best_new to combo chosen
        if proceed:
            combo_chosen.append(best_new)
        else:
            print('Search finally converged!')
    
    #Get out the optimal combo
    optimal_combo = combo_chosen[-1]
    #Format and save results
    results = pd.DataFrame(results)
    results.to_excel(path + 'forward_selector_results.xlsx')

    return(results, optimal_combo)


##


#########Example 
#signature = [40, 4, 25, 39, 18, 20, 10, 11, 17, 33]
#X_sig = X_all[:, signature]
#########

##


#sig_scorer(): given a gene signature, calculate each obs score.
def sig_scorer(X_sig):
    '''
    Given a signature (gene list) and an array of observations (with signature genes in columns)
    the function calculates each observation signature 'score' as the sum of each gene individual
    score, wich can be 1 or -1 (1 if >= median and -1 vice versa).
    '''

    #Select signature colums 
    medians = np.median(X_sig, axis=0)

    #Calculate scores for each observation
    scores = []
    for i in range(X_sig.shape[0]):
        s = []
        for j in range(X_sig.shape[1]):
            if X_sig[i, j] >= medians[j]:
                x = 1
            else:
                x = -1
            s.append(x)
        scores.append(sum(np.array(s)))
    
    return(np.array(scores))
    

##


#sig_scorer_cw(): continuos weighted alternative for signature scoring
def sig_scorer_cw(X_sig, y):
    '''
    Given a X_sig array, calculate a continuous weighted signature score.
    '''

    #Fit the model and retrieve coefficients
    model = CoxPHSurvivalAnalysis(n_iter = 1000000)
    model.fit(X_sig, y)
    weights = model.coef_
    #Calculate each obs weighted score
    scores = np.sum((weights / np.sum(weights)) * X_sig, axis = 1)

    return(scores)


##


#sig_scorer_cw(): continuos weighted alternative for signature scoring
def sig_scorer_dw(X_sig, y):
    '''
    Given a X_sig array, calculate a discrete weighted signature score.
    '''

    #Fit the model and retrieve coefficients
    model = CoxPHSurvivalAnalysis(n_iter = 1000000)
    model.fit(X_sig, y)
    weights = model.coef_

    #Calculate discrete indicators
    medians = np.median(X_sig, axis=0)
    X_d = np.ones((X_sig.shape[0], X_sig.shape[1]))
    for i in range(X_sig.shape[0]):
        for j in range(X_sig.shape[1]):
            if X_sig[i, j] >= medians[j]:
                X_d[i, j] = 1
            else:
                X_d[i, j] = -1
    #Calculate each obs weighted discrete score
    scores = np.sum((weights / np.sum(weights)) * X_d, axis = 1)

    return(scores)


#######Testing
#scores = sig_scorer(X[:, 100:105])
######  


##


#######Example
#path = '/Users/IEO5505/Desktop/ML_TDA/TCGA_data/plots_results/TCGA/enet_OS_TCGA_all/'
#######


#KM_plotter(): as the name suggests :)
def KM_Plotter(X_sig, y, kind, path):
    '''
    Given X filtered for genes of the signature obtained, an y target and a path, 
    plot the KM curve.
    '''
    list(y)
    #Calculate signature scores, and assign labels accordingly
    scores = sig_scorer(X_sig)
    label = np.array([ 'up' if x >= 0 else 'down' for x in scores ])

    #Estimate and compare survival curves of up and downs groups
    pval = compare_survival(y, label)[1]

    #Plot
    plt.figure(dpi=200.0)
    for value in np.unique(label): 
        mask = label == value
        time, p = kaplan_meier_estimator(y['delta'][mask], y['time'][mask])
        p = plt.step(time, p, where = "post", label = "%s (n = %d)" % (value, mask.sum()), 
        marker = '|', linewidth = 1.5)

    p = plt.ylabel("est. probability of survival $\hat{S}(t)$")
    p = plt.xlabel("time $t$ (days)")
    if kind == 'OS':
        p = plt.title('Kaplan-Meier curve: Overall survival', fontweight='bold')
    else:
        p = plt.title('Kaplan-Meier curve: RFS survival', fontweight='bold')
    p = plt.legend(loc = "best")
    p = plt.figtext(.68, .7, 'pvalue = %.3f' % pval)
    p = plt.savefig(path + 'KM_plot.png')
    p = plt.show()


#######Testing
#KM_Plotter(X_sig, y, path)
#######


##


########Example
#alphas = list(np.arange(0.0, 30, 0.1))
#ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#n_combos = 500
#
#alpha = alphas[20]
#ratio = ratios[-1]
########


##


#internal_evaluator_en(): same as before, but with en and a couple of input hyperparameters
def internal_evaluator_en(X, y, alpha, ratio, folds):
    '''
    Given an X and y array, some alpha and l1_ratio combo, the function evaluate by 5-fold CV the 
    predictive power of such a combination of hyperparameters, computing the C_index of the resulting 
    enCoxPH models.
    ''' 

    #Initialize C_score and hratio lists
    C = np.full(folds, 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over folds indeces. Fit on tr; score model on t; and store infos.
    fold_c = 0
    for tr, t in skf.split(np.zeros(len(y)), y['delta']): 
        m = CoxnetSurvivalAnalysis(alphas = [alpha], l1_ratio=ratio, max_iter=10000000)
        m.fit(X[tr, ], y[tr])
        C[fold_c] = m.score(X[t, ], y[t])
        fold_c += 1
    C_mean = np.mean(C)

    return(C_mean)


##


def rf_evaluator(X, y, n_est, min_samples_split, min_samples_leaf, max_f, folds):
    '''
    Given an X and y array, by 5-fold CV the predictive power of survival rf with
    a certain combo of hyperparameters.
    ''' 

    #Initialize C_score and hratio lists
    C = np.full(folds, 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over folds indeces. Fit on tr; score model on t; and store infos.
    fold_c = 0
    for tr, t in skf.split(np.zeros(len(y)), y['delta']): 
        m = RandomSurvivalForest(n_estimators=n_est,
                           min_samples_split=min_samples_split,
                           min_samples_leaf=min_samples_leaf,
                           max_features=max_f,
                           n_jobs=-1,
                           random_state=23)
        m.fit(X[tr, ], y[tr])
        C[fold_c] = m.score(X[t, ], y[t])
        fold_c += 1
    C_mean = np.mean(C)

    return(C_mean)


##


def gbCox_linear_evaluator_parallel(i, X, y, folds, combo):
    '''
    Given an X and y array, by 5-fold CV the predictive power of survival gradient 
    boosted linear Cox model with a certain combo of hyperparameters.
    ''' 

    #Initialize C_score and hratio lists
    C = np.full(folds, 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over folds indeces. Fit on tr; score model on t; and store infos.
    fold_c = 0
    for tr, t in skf.split(np.zeros(len(y)), y['delta']): 
        m = ComponentwiseGradientBoostingSurvivalAnalysis(learning_rate=combo[1], 
                                                        n_estimators=combo[0])
        m.fit(X[tr, ], y[tr])
        C[fold_c] = m.score(X[t, ], y[t])
        fold_c += 1
    C_mean = np.mean(C)

    return(i, C_mean)


##


def gbCox_trees_evaluator(X, y, folds, lr, n_est):
    '''
    Given an X and y array, by 5-fold CV the predictive power of survival gradient 
    boosted tree based Cox model with a certain combo of hyperparameters.
    ''' 

    #Initialize C_score and hratio lists
    C = np.full(folds, 1.111)

    #Define splitter for CV
    skf = StratifiedKFold(n_splits = folds, shuffle = True, random_state = 15)

    ##Loop over folds indeces. Fit on tr; score model on t; and store infos.
    fold_c = 0
    for tr, t in skf.split(np.zeros(len(y)), y['delta']): 
        m = GradientBoostingSurvivalAnalysis(learning_rate=lr, n_estimators=n_est, 
                                            max_depth=1)
        m.fit(X[tr, ], y[tr])
        C[fold_c] = m.score(X[t, ], y[t])
        fold_c += 1
    C_mean = np.mean(C)

    return(C_mean)


##


#parallel_GS_linear_GB(): Implements parallel GridSearch for tune the hyperparameters 
#of nonlinear models.
def parallel_GS_linear_GB(X, y, model, grid, n_combos, folds):
    '''
    Implements a parallel GridSearch for tune the hyperparameters of 
    ComponentwiseGradientBoostingSurvivalAnalysis.
    '''

    #Produce all combos
    all_combos = list(product( * [ list(grid[key]) for key in grid.keys() ]))
    
    #Sample n_combos
    random.seed(1234)
    combos = random.sample(all_combos, n_combos)

    #Perform grid search (parallelize the for loop)    
    with mp.Pool(8) as pool:
        result_objects = [ pool.apply_async(gbCox_linear_evaluator_parallel, args=(i, X, y, folds, combo)) for i, combo in enumerate(combos) ]
        C = [ r.get()[1] for r in result_objects ]
        pool.close()
        pool.join()

    #Choose the best combo
    best_C = max(C)
    best_combo = combos[C.index(best_C)]

    #Format results
    results = pd.DataFrame({'C_score' : C, 'hyperparameters' : combos})

    #Return
    return(best_C, best_combo, results)


##


############Example
#X = X[:, idx_features_to_use_all]
#X.shape
#
#model = ComponentwiseGradientBoostingSurvivalAnalysis(subsample=1)
#
#grid = {
#    'n_learners' : list(np.arange(50, 1000, 20)),
#    'learning_rate' : [0.0001, 0.005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]
#    }
#
#n_combos = 10
############

############Testing 
#a = timer()
#best_C_p, best_combo_p, results_p = parallel_GS_linear_GB(X, y, model, grid, n_combos, folds)
#b = timer()
#time_parallel = b - a 
###########


##


#########Example 
signature = ['LDHA', 'EIF2S2', 'KRT7', 'CSDE1', 'SEC14L1', 'MIR4435-2HG', 
    'DDIT4', 'HMGB3', 'FN1', 'DGKD', 'TTC3', 'SQSTM1', 
    'HES1', 'ZNF217', 'DDX6', 'FTL', 'SLC7A11']
task = 'OS'
score_type = 'dd'
subsetted = False
path = '/Users/IEO5505/Desktop/'
########


##


#sig_evaluator(): implements multivariate coxphmodel evaluation of a signature score 
#on the training set from which it has been chosen
def sig_evaluator(df, y, signature, task, score_type, path, subsetted = False):
    '''
    Given a database dataframe, a signature to test, a survival task and a score_type, 
    the function fits a multivariate CoxRegression model with a signature score and the clinical 
    features taken into account in the study, to test for the impact of the brand new molecular
    signature on patients survival, given other known prognostic factors.
    '''

    #Prepare the dataframe and produce the gene expression array for sig score calculation
    
    ##Decide task
    if task == 'OS':
        surv_columns = ['OS_status', 'OS_time'] 
    else:
        surv_columns = ['Prog_free_status', 'Prog_free_time'] 
    ##Select columns
    if subsetted == True:
        clinical = ['age_at_diagnosis', 'stage']
    else:
        clinical = ['age_at_diagnosis', 'stage', 'subtype']
    cols = signature + clinical + surv_columns
    df_m = df.loc[:, cols]
    ##Rename and reformat columns
    if task == 'OS': 
        df_m = df_m.rename(columns={"OS_status": "Status", "OS_time": "Time"})
    else:
        df_m = df_m.rename(columns={"Prog_free_status": "Status", "Prog_free_time": "Time"})
    df_m['Status'] = df_m['Status'].astype(int)
    for col in [ x for x in clinical if x != 'age_at_diagnosis' ]: 
        df_m[col] = df_m[col].astype('category')
    df_m = OneHotEncoder().fit_transform(df_m)
    #Produce X_sig, adn drop signature columns
    X_sig = df_m.loc[:, signature].values
    df_m = df_m.drop(signature, axis = 1)

    #Add a sig_score column
    if score_type == 'discrete':
        df_m['sig_score'] = sig_scorer(X_sig)
    elif score_type == 'discrete_weighted':
        df_m['sig_score'] = sig_scorer_dw(X_sig, y)
    else:
        df_m['sig_score'] = sig_scorer_cw(X_sig, y)

    #Fit the lifelines CoxPHmodel
    cph = CoxPHFitter()
    cph.fit(df_m, 'Time', event_col='Status')

    #Make the HRs plot and save the model summary:

    ##Take the HRs and the CI of each feature, and of the signature score
    HR = cph.hazard_ratios_
    CI = cph.confidence_intervals_
    pval = cph.summary.loc['sig_score', 'p']
    hr = HR['sig_score']
    l_ci = np.exp(CI).loc['sig_score', '95% lower-bound']
    u_ci = np.exp(CI).loc['sig_score', '95% upper-bound']
    #Make the plot
    plt.figure(dpi=220, figsize=(5, 4))
    p = cph.plot(hazard_ratios=True)
    p = plt.title('HR: OS on all METABRIC', fontweight='bold')
    p = plt.figtext(.65, .25, 'HR = %.3f (CI: %.3f-%.3f) \np = %.3f' % (hr, l_ci, u_ci, pval), fontsize = 'xx-small')
    p = plt.subplots_adjust(left=0.4, bottom=0.2)
    p = plt.savefig(path + 'HR_plot.png')
    p = plt.show()

    #Save and print the summary
    cph.summary.to_excel(path + 'Multivariate_model_summary.xlsx')
    cph.print_summary()


#####Testing
#sig_evaluator(df, y, signature, task, score_type, path, subsetted = False)
#####


##
