"""
Integrated analysis utils.
"""

import numpy as np
import pandas as pd


##


# Utils
def filter_loops_with_enhancers(df, df_enhancers):
    """
    Test if as set of enhancers coincide with some HiChip loop (i.e., HiChip bin interaction pair).
    """
    L = []
    for i in range(df_enhancers.shape[0]):
        d = df_enhancers.iloc[i].to_dict()
        enh_chr = d['chr']
        enh_pos = d['start']
        cluster = d['cluster']
        test1 = (enh_pos >= df['start1']) & (enh_pos <= df['end1']) & (enh_chr == df['seqnames1'])
        test2 = (enh_pos >= df['start2']) & (enh_pos <= df['end2']) & (enh_chr == df['seqnames2'])
        if test1.any():
            isin = 1
            idx = np.where(test1)[0]
        elif test2.any():
            isin = 2
            idx = np.where(test2)[0]
        start1, end1 = df.iloc[idx,1:3].values[0]
        start2, end2 = df.iloc[idx,4:6].values[0]
        l = [enh_chr, start1, end1, enh_chr, start2, end2, enh_pos, cluster, isin]
        L.append(l)
    
    # DataFrame
    df_filtered = pd.DataFrame(
        L,
        columns=[
            'seqnames1', 'start1', 'end1', 
            'seqnames2', 'start2', 'end2', 
            'enhancer_summit', 'cluster', 'isin_bin'
        ]
    )
    
    return df_filtered
        

##


def filter_loops_with_tss(df, df_tss):
    """
    Filter bedpe file with tss table.
    """
    L = []
    for i in range(df.shape[0]):
        d = df.iloc[i].to_dict()
        chrom = d['seqnames1']
        isin_bin = d['isin_bin']
        start1 = d['start1']
        end1 = d['end1']
        start2 = d['start2']
        end2 = d['end2']
        if isin_bin == 0:
            start = start1
            end = end1
        elif isin_bin == 1:
            start = start2
            end = end2
        test = (df_tss['chr'] == chrom) & (df_tss['tss'] >= start) & (df_tss['tss'] <= end)
        if test.any():
            gene_record = df_tss.loc[test]
            d['tss'] = gene_record['tss'].iloc[0]
            d['gene'] = gene_record['gene_name'].iloc[0]
            d['gene_bin'] = int(not isin_bin)
            L.append(list(d.values()))

    # DataFrame
    df_filtered = pd.DataFrame(
        L, 
        columns=[
            'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', 
            'enhancer_summit', 'enhancer_cluster', 'enh_bin', 'tss', 'gene', 'tss_bin'
        ]
    )

    return df_filtered


##
