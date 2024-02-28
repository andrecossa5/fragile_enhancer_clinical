#!/usr/bin/python

"""
Script for survival analysis of DE genes.
"""

import os
from utils.survival import *
from plotting_utils._plotting_base import *

##


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'survival')


##


def main():

    # Read formatted METABRIC data and DEG list 
    df = pd.read_csv(os.path.join(path_data, 'METABRIC', 'formatted_data.csv'), index_col=0)
    DEGs = pd.read_csv(os.path.join(path_data, 'expression', 'GRHL2_KD_DEG.csv'), index_col=0)
    DEGs = DEGs.loc[DEGs['FDR']<.1,['logFC', 'FDR']].sort_values('FDR').copy() # Remove FDR >= .1
    clinical_features = df.columns[:10].tolist()
    expression_features =  df.columns[df.columns.isin(DEGs.index)].tolist()

    ##

    # df = df.loc[df['OS_STATUS'] & df['RFS_STATUS']]
    # np.mean(df['OS_time']-df['RFS_time'])

    #


    # Perform RFS survival analysis
    status_cov = 'RFS_STATUS'
    time_cov  = 'RFS_time'
    confounders = ['AGE_AT_DIAGNOSIS', 'TUMOR_STAGE', 'subtype']

    # Univariate
    L = [
        univariate_Cox(
        df, molecular_feature=x, status_cov=status_cov, time_cov=time_cov) \
        for x in expression_features
    ]
    df_univariate = pd.concat(L, axis=1).T.set_index('molecular_feature')

    # Multivariate
    L = [
        multivariate_Cox(
        df, molecular_feature=x, status_cov=status_cov, time_cov=time_cov, confounders=confounders) \
        for x in expression_features
    ]
    df_multivariate = pd.concat(L, axis=1).T.set_index('molecular_feature')

    # Save
    df_results = pd.concat([df_univariate, df_multivariate])
    df_results.to_csv(os.path.join(path_results, 'survival_df.csv'))


    ##


    # Load and visualize with KMs.
    df_results = pd.read_csv(os.path.join(path_results, 'survival_df.csv'), index_col=0)
    df_results = df_results.join(DEGs)


    ##

    # Distribution time
    fig, axs = plt.subplots(1,2,figsize=(12,5))

    sns.kdeplot(df, x='OS_time', ax=axs[0], color='k', fill=True)
    format_ax(axs[0], xlabel=f'Overrall Survival (OS) t (days)', ylabel='Density',
            reduce_spines=True, yticks='')
    t = f'Cases: {df["OS_STATUS"].sum()} dead, {(~df["OS_STATUS"]).sum()} alive'
    mean_OS_t = f'Median OS time: {df["OS_time"].median():.2f} days'
    axs[0].text(.5, .9, t, transform=axs[0].transAxes)
    axs[0].text(.5, .85, mean_OS_t, transform=axs[0].transAxes)

    sns.kdeplot(df, x='RFS_time', ax=axs[1], color='k', fill=True)
    format_ax(axs[1], xlabel=f'Overrall Survival (RFS) t (days)', ylabel='Density', 
            reduce_spines=True, yticks='')
    t = f'Cases: {df["RFS_STATUS"].sum()} relapsed, {(~df["RFS_STATUS"]).sum()} not relapsed'
    mean_OS_t = f'Median RFS time: {df["RFS_time"].median():.2f} days'
    axs[1].text(.5, .9, t, transform=axs[1].transAxes)
    axs[1].text(.5, .85, mean_OS_t, transform=axs[1].transAxes)

    fig.suptitle('METABRIC cohort and survival problems')
    fig.subplots_adjust(left=.1, right=.9, bottom=.1, top=.9)
    fig.savefig(os.path.join(path_results, 'summary_problems.png'), dpi=300)


    ##


    # Univariate
    df_ = df_results.query('model_type == "univariate"')

    fig, axs = plt.subplots(1,3,figsize=(13,5))

    stem_plot(df_.sort_values('HR_mean', ascending=False).head(20), x='HR_mean', ax=axs[0])
    format_ax(axs[0], xlabel='Mean Hazard Ratio', reduce_spines=True)
    stem_plot(df_.sort_values('C_mean', ascending=False).head(20), x='C_mean', ax=axs[1])
    format_ax(axs[1], xlabel='Mean Harrel\'s C-index', reduce_spines=True)
    stem_plot(df_.sort_values('logFC', ascending=False).head(20), x='logFC', ax=axs[2])
    format_ax(axs[2], xlabel='Differential Expression logFC', reduce_spines=True)
    fig.suptitle(f'Relapse-Free Survival (RFS) top predictors vs DE top genes')

    fig.subplots_adjust(left=.1, right=.9, bottom=.1, top=.9, wspace=.4)
    fig.savefig(os.path.join(path_results, 'summary_univariate.png'), dpi=300)


    ##


    # Multivariate
    df_ = df_results.query('model_type == "multivariate"')

    fig, axs = plt.subplots(1,3,figsize=(13,5))

    stem_plot(df_.sort_values('HR_mean', ascending=False).head(20), x='HR_mean', ax=axs[0])
    format_ax(axs[0], xlabel='Mean Hazard Ratio', reduce_spines=True)
    stem_plot(df_.sort_values('C_mean', ascending=False).head(20), x='C_mean', ax=axs[1])
    format_ax(axs[1], xlabel='Mean Harrel\'s C-index', reduce_spines=True)
    stem_plot(df_.sort_values('logFC', ascending=False).head(20), x='logFC', ax=axs[2])
    format_ax(axs[2], xlabel='Differential Expression logFC', reduce_spines=True)
    fig.suptitle(f'Relapse-Free Survival (RFS) top predictors vs DE top genes')

    fig.subplots_adjust(left=.1, right=.9, bottom=.1, top=.9, wspace=.4)
    fig.savefig(os.path.join(path_results, 'summary_multivariate.png'), dpi=300)


    ##


    # DEGs
    up = df_results['logFC'].sort_values().tail(50).index.tolist()
    down = df_results['logFC'].sort_values().head(50).index.tolist()

    df_results['gene_status'] = np.select(
        [df_results.index.isin(up), df_results.index.isin(down)], 
        ['up', 'down'],
        default='other'
    )
    df_results['gene_status'] = pd.Categorical(
        df_results['gene_status'], categories=['up', 'down', 'other']
    )

    fig, axs = plt.subplots(1,2,figsize=(10,5))
    pairs = [['up', 'down'], ['up', 'other'], ['down', 'other']]
    box(
        df_results.query('model_type == "multivariate"'), 
        'gene_status', 'HR_mean', c='#915555', ax=axs[0], with_stats=True, pairs=pairs
    )
    format_ax(
        axs[0], title='RFS Hazard ratio', ylabel='Mean HR (5-fold CV)', 
        reduce_spines=True)
    box(
        df_results.query('model_type == "multivariate"'), 
        'gene_status', 'C_mean', c='#915555', ax=axs[1], with_stats=True, pairs=pairs
    )
    format_ax(
        axs[1], title='RFS C-index', ylabel='Mean C-index (5-fold CV)', 
        reduce_spines=True)

    fig.tight_layout()
    fig.savefig(os.path.join(path_results, 'top_and_boottom_DEGs.png'), dpi=300)


    ## 


    # Confounders relevance

    # Top genes
    df_top = df_results.loc[df_results.index.isin(up) | df_results.index.isin(down)]
    (
        df_top.query('HR_mean>=1. and C_mean>.5')
        .to_csv(os.path.join(path_results, 'interesting_genes.csv'))
    )



    t = df_results.query('model_type == "univariate"')['HR_mean'].sort_values().tail(50).index

    set(t) & set(down)


    ##


# Run
if __name__ == '__main__':
    main()
