# Loops

# Code
library(tidyverse)
library(HiCDCPlus)

# Paths
path_pairs <- '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/Hi_Chip/hic_pro/SCR_pool_corr.txt'
outdir <- '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/Hi_Chip/hic_pro'

# Generate genomic features
construct_features(
  output_path=paste0(outdir,"/hg19_5kb_GATC"),
  gen="Hsapiens", gen_ver="hg19",
  sig="GATC", bin_type="Bins-uniform",
  binsize=5000
)

# Create and fill GI list
gi_list <- generate_bintolen_gi_list(bintolen_path=paste0(outdir, '/hg19_5kb_GATC_bintolen.txt.gz'))
gi_list <- add_hicpro_allvalidpairs_counts(gi_list, path_pairs)
gi_list <- expand_1D_features(gi_list)

# HiCDCPlus
set.seed(1234)
gi_list <- HiCDCPlus_parallel(gi_list, ncore=8)
head(gi_list)
gi_list_write(gi_list, fname=paste0(outdir,'/loop_calling.txt.gz'))


##
