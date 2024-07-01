
library(tidyverse)
library(ggpubr)

SEED <- 4321
set.seed(SEED)

path_output_temp <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/temp/")

SOURCES <- c("icgc", "hart")
MARKERS <- c("CtIP", "GRHL")
WIN <- 50
Mb <- 1000000

high_low <- T


##


# Define groups of clusters
cluster_groups <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
cluster_groups$CtIP <- list(
  "high" = c("CtIP_cluster_2_Enh", "CtIP_cluster_3_Enh", "CtIP_cluster_5_Enh"), 
  "low" = c("CtIP_cluster_6.2_Enh", "CtIP_cluster_6.3_Enh", "CtIP_cluster_6.1_Enh", "CtIP_cluster_6.0")
)
cluster_groups$GRHL <- list(
  "high" = c("GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh"), 
  "low" = c("GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.0")
)

# Read FC distributions from ICGC and Hartwig
fc_all_samples_markers_source <- list()
for(source in SOURCES){
  fc_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.all_enhancers.", source, ".rds")))   
  fc_ran_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.all_random.", source, ".rds")))  
  fc_x_group_all_samples_markers_source[[source]] <- readRDS(file = fs::path(path_output_temp, paste0("fc_dist.grouped_enhancers.", source, ".rds")))   
}


##


## PLOT

for(marker in MARKERS){
  
}