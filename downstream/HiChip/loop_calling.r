# Loop calling

# Code
library(tidyverse)
library(HiCDCPlus)

# Paths
# path_pairs <- '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/Hi_Chip/hic_pro/SCR_pool.txt'
# outdir <- '/Users/IEO5505/Desktop/fragile_enhancer_clinical/data/Hi_Chip/hic_pro'
# res <- as.numeric("50000")

# Parse args
args <- commandArgs(TRUE)
path_pairs <- args[1] 
outdir <- args[2]
res <- as.numeric(args[3])
sample_name <- args[4]

# Generate genomic features
construct_features(
  output_path=paste0(outdir, '/', round(res/1000), '_kb'),
  gen="Hsapiens", gen_ver="hg19",
  sig="GATC", bin_type="Bins-uniform",
  binsize=res
)

# Create and fill GI list
gi_list <- generate_bintolen_gi_list(
  bintolen_path=paste0(outdir, '/', round(res/1000), '_kb_bintolen.txt.gz')
)
gi_list <- add_hicpro_allvalidpairs_counts(gi_list, path_pairs)
gi_list <- expand_1D_features(gi_list)

# HiCDCPlus
set.seed(1234)
gi_list <- HiCDCPlus_parallel(gi_list, ncore=2)
gi_list_write(
  gi_list, 
  fname=paste0(outdir, '/', sample_name, '_loops_', round(res/1000), '.txt.gz')
)


#
