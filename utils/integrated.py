"""
Integrated analysis utils.
"""

import numpy as np
import pandas as pd


# Utils
def filter_enhancers(enh, df):
    """
    Test if as set of enhancers coincide with some HiChip loop (i.e., HiChip bin interaction pair).
    """
    L = []
    for i in range(enh.shape[0]):
        d = enh.iloc[i].to_dict()
        enh_chr = d['chrI']
        enh_pos = d['startI']
        cluster = d['cluster']
        testI = (enh_pos >= df['startI']) & (enh_pos <= df['endI']) & (enh_chr == df['chrI'])
        testJ = (enh_pos >= df['startJ']) & (enh_pos <= df['endJ']) & (enh_chr == df['chrJ'])
        if testI.any():
            isin = 0
            idx = np.where(testI)[0]
        elif testJ.any():
            isin = 1
            idx = np.where(testJ)[0]
        start_i, end_i = df.iloc[idx,1:3].values[0]
        start_j, end_j = df.iloc[idx,4:6].values[0]
        l = [enh_chr, start_i, end_i, enh_chr, start_j, end_j, enh_pos, cluster, isin]
        L.append(l)
    
    # DataFrame
    df_coords = pd.DataFrame(L, columns=['chr', 'startI', 'endI', 'chr', 'startJ', 'endJ', 'enhancer_summit', 'cluster', 'isin_bin'])
    
    return df_coords
        

##


def filter_tss(df_tss, df):
    """
    Filter bedpe file with tss table
    """
    L = []
    for i in range(df.shape[0]):
        d = df.iloc[i].to_dict()
        chrom = d['chr']
        isin_bin = d['isin_bin']
        start_i = d['startI']
        end_i = d['endI']
        start_j = d['startJ']
        end_j = d['endJ']
        if isin_bin == 0:
            start = start_i
            end = end_i
        elif isin_bin == 1:
            start = start_j
            end = end_j
        test = (df_tss['chr'] == chrom) & (df_tss['tss'] >= start) & (df_tss['tss'] <= end)
        if test.any():
            gene_record = df_tss.loc[test]
            d['tss'] = gene_record['tss'].iloc[0]
            d['gene'] = gene_record['gene_name'].iloc[0]
            d['gene_bin'] = int(not isin_bin)
            L.append(list(d.values()))

    # DataFrame
    df_coords = pd.DataFrame(
        L, columns=[
            'chr', 'startI', 'endI', 'startJ', 'endJ', 'enhancer_summit', 
            'enhancer_cluster', 'enh_bin', 'tss', 'gene', 'tss_bin'
        ]
    )

    return df_coords


##
