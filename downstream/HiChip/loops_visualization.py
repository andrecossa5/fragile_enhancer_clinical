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
path_cool = os.path.join(path_main, 'data', 'functional_genomics', 'HiChip', 'mcool')
path_results = os.path.join(path_main, 'results', 'integrated')


##


# Read .mcool and regions

# SCR
sample = 'SCR'
# cooler.fileops.list_coolers(os.path.join(path_cool, f'{sample}.mcool'))
C_SCR = { 
    '2000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/2000')), 
    '4000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/4000')), 
    '8000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/8000')), 
}

# KD
sample = 'KD'
C_KD = { 
    '2000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/2000')), 
    '4000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/4000')), 
    '8000' : cooler.Cooler(os.path.join(path_cool, f'{sample}.mcool::resolutions/8000')), 
}


# Read loops
df_loops = pd.read_csv(os.path.join(path_results, 'prova_loops.csv'), index_col=0)


##


# Plot loops

# 1
title = f'CtIP enhancer cluster 2 (summit chr8:{102449518}) vs GRHL2 TSS (chr8:{102504667})'
chrom = '8'
start = 102448000
end = 102508000 
res = '4000'
span = 18000
vmax = .1
balance = True
lognorm = True
norm = 'ICE'

fig, axs = plt.subplots(1,2,figsize=(10,5))
plot_region(C_SCR[res], chrom, start-span, end+span, balance=balance, lognorm=lognorm, vmax=vmax, ax=axs[0]) 
plot_region(C_KD[res], chrom, start-span, end+span, balance=balance, lognorm=lognorm, vmax=vmax, ax=axs[1])
fig.tight_layout()
fig.suptitle(f'{title}, {res} bp')
fig.savefig(os.path.join(path_results, f'{title}_{res}_kb_{norm}.png'), dpi=300)


##


# 2
title = f'CtIP enhancer cluster 2 (summit chr1:{201277106}) vs LAD1 TSS (chr1:{201368452})'
chrom = '1'
start = 201276000
end = 201372000 
res = '4000'
span = 30000
vmax = .1
balance = True
lognorm = True
norm = 'ICE'

fig, axs = plt.subplots(1,2,figsize=(10,5))
plot_region(C_SCR[res], chrom, start-span, end+span, balance=balance, lognorm=lognorm, vmax=vmax, ax=axs[0]) 
plot_region(C_KD[res], chrom, start-span, end+span, balance=balance, lognorm=lognorm, vmax=vmax, ax=axs[1])
fig.tight_layout()
fig.suptitle(f'{title}, {res} bp')
fig.savefig(os.path.join(path_results, f'{title}_{res}_kb_{norm}.png'), dpi=300)


##