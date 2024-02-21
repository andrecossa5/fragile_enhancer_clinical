"""
Exploratory maps. Hi-chip.
"""

import os
import cooler
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from cooltools.lib import plotting
matplotlib.use('macOSX')
matplotlib.rcParams["axes.formatter.useoffset"] = False
matplotlib.rcParams["axes.formatter.use_mathtext"] = False


## 


# Utils
def ticks_formatter(ax, x=True, y=True, chrom='#', res=1000, rotate=True):
    """
    Format axes ticks and labels.
    """

    from matplotlib.ticker import EngFormatter
    bp_formatter = EngFormatter('b')

    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    else:
        ax.set_yticks([])
        ax.yaxis.set_visible(False)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    else:
        ax.set_xticks([])
        ax.xaxis.set_visible(False)
    if rotate:
        ax.tick_params(axis='x',rotation=45)

    ax.set_xlabel(f'chr{chrom} ({res} kb binsize)')


##


def plot_region(c, chrom, start, end, lognorm=False,
    title=None, cmap='fall', vmax=10, balance=False, ax=None,
    show_x=True, show_y=True, rotate=True
    ):
    """
    Plot a small region inside a chromosome.
    """

    region = (chrom, start, end)
    M = c.matrix(balance=balance).fetch(region)

    from matplotlib.colors import LogNorm
    norm = LogNorm(vmax=vmax)
    norm = norm if (not balance and lognorm) else None
    vmax = vmax if norm is None else None

    im = ax.matshow(
        M, extent=(start, end, end, start),
        norm=norm, vmax=vmax, cmap=cmap
    )
    ax.set(title=title)
    value_type = 'balanced' if balance else 'raw'
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(value_type)
    ticks_formatter(ax, chrom=chrom, res=c.binsize, 
                    x=show_x, y=show_y, rotate=rotate)

    return ax


##


# Paths
path_main = '/Users/IEO5505/Desktop/fragile_enhancer_clinical'
path_data = os.path.join(path_main, 'data', 'Hi_Chip', 'cool')
path_results = os.path.join(path_main, 'results', 'Hi_Chip')


##


# Read .mcool and regions
sample = 'KD_pool_120'

cooler.fileops.list_coolers(os.path.join(path_data, f'{sample}.mcool'))
C = { 
    '2000' : cooler.Cooler(os.path.join(path_data, f'{sample}.mcool::resolutions/2000')), 
    '8000' : cooler.Cooler(os.path.join(path_data, f'{sample}.mcool::resolutions/8000')), 
    '16000' : cooler.Cooler(os.path.join(path_data, f'{sample}.mcool::resolutions/16000')), 
    '32000' : cooler.Cooler(os.path.join(path_data, f'{sample}.mcool::resolutions/32000')), 
}

# Regions
regions = {
    'EMP2' : {
        'chrom':'16', 
        'start' : 10600000, # -500000, 
        'end' : 10700000, # +500000 
    },
}


##


# Here we go

# 1. Regions at with different methods, highest resolution
for gene in regions:

    chrom = regions[gene]['chrom']
    start = regions[gene]['start']
    end = regions[gene]['end']

    res = '2000'
    vmax_raw = 3
    vmax_balanced = .01

    fig, axs = plt.subplots(2,2,figsize=(10,7))

    plot_region(C[res], chrom, start, end, title='Raw, linear', balance=False, lognorm=False, 
            vmax=vmax_raw, ax=axs[0,0], show_x=False)
    plot_region(C[res], chrom, start, end, title='Raw, lognorm', balance=False, lognorm=True, 
            vmax=vmax_raw, ax=axs[0,1], show_x=False, show_y=False)
    plot_region(C[res], chrom, start, end, title='Balanced, linear', balance=True, lognorm=False,
            vmax=vmax_balanced, ax=axs[1,0])
    plot_region(C[res], chrom, start, end, title='Balanced, lognorm', balance=True, lognorm=True, 
            vmax=vmax_balanced, ax=axs[1,1], show_y=False)

    fig.subplots_adjust(top=.9, bottom=.15, left=.2, right=.8, wspace=.3, hspace=.2)
    fig.suptitle(f'{gene} locus, binsize {res} kb')
    
    fig.savefig(os.path.join(path_results, f'{sample}_{gene}_{res}.png'), dpi=300)


##


# 2. Regions 4 resolutions
for gene in regions:

    chrom = regions[gene]['chrom']
    start = regions[gene]['start']
    end = regions[gene]['end']

    vmax_raw = 3
    vmax_balanced = .01

    fig, axs = plt.subplots(2,4,figsize=(16,7))

    for i, res in enumerate(C):
        plot_region(
            C[res], chrom, start, end, title=f'Raw, {res} kb', balance=False, 
            lognorm=False, vmax=vmax_raw*(i+1), ax=axs[0,i], 
            show_x=False, 
            show_y=False if i != 0 else True
        )

    for i, res in enumerate(C):
        plot_region(
            C[res], chrom, start, end, title=f'Balanced, {res} kb', balance=True, 
            lognorm=False, vmax=vmax_balanced*(i+1), ax=axs[1,i], 
            show_x=True, 
            show_y=False if i != 0 else True
        )

    fig.subplots_adjust(top=.9, bottom=.15, left=.1, right=.9, wspace=.4, hspace=.1)
    fig.suptitle(f'{gene} locus, different resolution')
    fig.savefig(os.path.join(path_results, f'{sample}_{gene}.png'), dpi=300)


##