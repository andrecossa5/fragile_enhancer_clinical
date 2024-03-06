"""
Exploratory maps. Hi-chip.
"""

import matplotlib
import matplotlib.pyplot as plt
from cooltools.lib import plotting


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

    ax.set_xlabel(f'chr{chrom} ({int(res/1000)} kb)')


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