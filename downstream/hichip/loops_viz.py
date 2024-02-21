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
path_cool = os.path.join(path_main, 'data', 'Hi_Chip', 'cool')
path_loops = os.path.join(path_main, 'data', 'Hi_Chip', 'loops')
path_results = os.path.join(path_main, 'results', 'Hi_Chip')


##


# Read .mcool and regions
sample = 'KD_pool_120'
cooler.fileops.list_coolers(os.path.join(path_cool, f'{sample}.mcool'))
C = { 
    '2000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/2000')), 
    '8000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/8000')), 
    '16000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/16000')), 
    '32000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/32000')), 
}

# Read loops
df_loops = pd.read_csv(os.path.join(path_results, 'SCR_CtIP_enh_loops.csv'), index_col=0)


##


# Plot loop
# chr20_352000_354000_chr_20_388000_390000
chrom = '20'
start = 350000
end = 392000
res = '2000'
vmax_raw = 30
vmax_balanced = .02

fig, ax = plt.subplots(figsize=(7,7))
plot_region(C[res], chrom, start, end, balance=False, lognorm=False, vmax=vmax_raw, ax=ax)
fig.tight_layout()
plt.show()
# fig.subplots_adjust(top=.9, bottom=.15, left=.2, right=.8, wspace=.3, hspace=.2)
# fig.suptitle(f'{gene} locus, binsize {res} kb')
# fig.savefig(os.path.join(path_results, f'{sample}_{gene}_{res}.png'), dpi=300)


##
