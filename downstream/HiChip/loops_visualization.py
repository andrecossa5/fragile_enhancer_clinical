"""
Visualize loops.
"""

import os
import cooler
import numpy as np
import pandas as pd
from utils.hichip_visualization import *
matplotlib.use('macOSX')
matplotlib.rcParams["axes.formatter.useoffset"] = False
matplotlib.rcParams["axes.formatter.use_mathtext"] = False


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
