
### CtIP/GRHL2 enhancers overlap ###

library(tidyverse)
library(GenomicRanges)

SEED <- 4321
set.seed(SEED)

path_main <- "/Users/ieo6983/Desktop/fragile_enhancer_clinical/"
path_enhancers <- fs::path(path_main, "data/functional_genomics/others/")
path_results <- fs::path(path_main, "results/Chip/")


##

MARKERS <- c("CtIP", "GRHL2")

# Read enhancers 
enh_all <- vector(mode = "list", length = 2)
columns_names <- c("chrom", "start", "end", "cluster")
enh_all <- lapply(list("Cluster_CtIP_Enh_All.txt", "Cluster_GRHL_Enh_All.txt"), FUN = function(file_name){fs::path(path_enhancers, file_name)}) %>%
  lapply(., function(x){read_tsv(x, col_names = columns_names, comment = "#")})
names(enh_all) <- c("CtIP", "GRHL2")

# add summit & name
enh_all <- lapply(enh_all, function(df){
  df$summit <- df$end
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})


##


## Add comparable and sensible cluster information 
# Creating clusters conversion table
df_comp <- data.frame("cluster_orig" = c(unique(enh_all$CtIP$cluster), unique(enh_all$GRHL2$cluster)), 
                      "cluster_new" = "tbd", 
                      "group" = "tbd",
                      "marker" = c(
                        rep("CtIP", length(unique(enh_all$CtIP$cluster))), 
                        rep("GRHL2", length(unique(enh_all$GRHL2$cluster)))
                      ))

# CtIP clusters
df_comp[df_comp$cluster_orig == "CtIP_cluster_2_Enh", ]$cluster_new <- "1"
df_comp[df_comp$cluster_orig == "CtIP_cluster_3_Enh", ]$cluster_new <- "2"
df_comp[df_comp$cluster_orig == "CtIP_cluster_5_Enh", ]$cluster_new <- "3"
df_comp[df_comp$cluster_orig == "CtIP_cluster_6.2_Enh", ]$cluster_new <- "4"
df_comp[df_comp$cluster_orig == "CtIP_cluster_6.3_Enh", ]$cluster_new <- "5"
df_comp[df_comp$cluster_orig == "CtIP_cluster_6.1_Enh", ]$cluster_new <- "6"
df_comp[df_comp$cluster_orig == "CtIP_cluster_6.0", ]$cluster_new <- "7"

# GRHL2 clusters 
df_comp[df_comp$cluster_orig == "GRHL_cluster_1_Enh", ]$cluster_new <- "1"
df_comp[df_comp$cluster_orig == "GRHL_cluster_2_Enh", ]$cluster_new <- "2"
df_comp[df_comp$cluster_orig == "GRHL_cluster_4_Enh", ]$cluster_new <- "3"
df_comp[df_comp$cluster_orig == "GRHL_cluster_5.1_Enh", ]$cluster_new <- "4"
df_comp[df_comp$cluster_orig == "GRHL_cluster_5.2_Enh", ]$cluster_new <- "5"
df_comp[df_comp$cluster_orig == "GRHL_cluster_5.3_Enh", ]$cluster_new <- "6"
df_comp[df_comp$cluster_orig == "GRHL_cluster_5.0", ]$cluster_new <- "7"

# High and Low clusters
df_comp[df_comp$cluster_new %in% c("1", "2", "3"), ]$group <- "high"
df_comp[df_comp$cluster_new %in% c("4", "5", "6", "7"), ]$group <- "low"

#df_comp %>% write_tsv(., fs::path(path_results, "Clusters_conversion.tsv"))

enh_all <- enh_all %>% lapply(., function(x){
  left_join(x, df_comp[, c("cluster_orig", "cluster_new", "group")], 
            by = c("cluster" = "cluster_orig"))
})


##


## Overlap among extended enhancers
WIN <- 50 

enh_all_ext <- lapply(enh_all, function(x){
  print(paste0("Extending enhancer by +/- ", WIN, "bp"))
  x$start <- x$summit - WIN
  x$end <- x$summit + WIN
  x <- makeGRangesFromDataFrame(x, keep.extra.columns = T)
  return(x)
}) 

for(clust in unique(enh_all_ext$CtIP$cluster_new)){
  print(paste0("Checking overlaps in cluster: ", clust))
  
  ctip_clust <- enh_all_ext$CtIP[enh_all_ext$CtIP$cluster_new == clust]
  grhl2_clust <- enh_all_ext$GRHL2[enh_all_ext$GRHL2$cluster_new == clust]
  
  ovrlps_ctip <- findOverlaps(ctip_clust, grhl2_clust, select = "arbitrary")
  ovrlps_grhl2 <- findOverlaps(grhl2_clust, ctip_clust, select = "arbitrary")
  
  # How many CtIP enhancers overlap with a GRHL2 enhancer of the same cluster?
  n_overlaps_ctip <- ctip_clust[!is.na(ovrlps_ctip)] %>% length
  n_ctip_tot <- length(ctip_clust)
  overlaps_ctip_perc <- (n_overlaps_ctip / n_ctip_tot)*100
  
  # Viceversa
  n_overlaps_grhl2 <- grhl2_clust[!is.na(ovrlps_grhl2)] %>% length
  n_grhl_tot <- length(grhl2_clust)
  overlaps_grhl2_perc <- (n_overlaps_grhl2 / n_grhl_tot)*100
  
  print(paste0("Percentage of CtIP enhancers with an overlapping GRHL2 enhancer: ", round(overlaps_ctip_perc), "%"))
  print(paste0("Percentage of GRHL2 enhancers with an overlapping CtIP enhancer: ", round(overlaps_grhl2_perc), "%"))
  
}

# "high" vs. "low" enhancers
for(clust in unique(enh_all_ext$CtIP$group)){
  print(paste0("Checking overlaps in cluster: ", clust))
  
  ctip_clust <- enh_all_ext$CtIP[enh_all_ext$CtIP$group == clust]
  grhl2_clust <- enh_all_ext$GRHL2[enh_all_ext$GRHL2$group == clust]
  
  ovrlps_ctip <- findOverlaps(ctip_clust, grhl2_clust, select = "arbitrary")
  ovrlps_grhl2 <- findOverlaps(grhl2_clust, ctip_clust, select = "arbitrary")
  
  # How many CtIP enhancers overlap with a GRHL2 enhancer of the same cluster?
  n_overlaps_ctip <- ctip_clust[!is.na(ovrlps_ctip)] %>% length
  n_ctip_tot <- length(ctip_clust)
  overlaps_ctip_perc <- (n_overlaps_ctip / n_ctip_tot)*100
  
  # Viceversa
  n_overlaps_grhl2 <- grhl2_clust[!is.na(ovrlps_grhl2)] %>% length
  n_grhl_tot <- length(grhl2_clust)
  overlaps_grhl2_perc <- (n_overlaps_grhl2 / n_grhl_tot)*100
  
  print(paste0("Percentage of CtIP enhancers with an overlapping GRHL2 enhancer: ", round(overlaps_ctip_perc), "%"))
  print(paste0("Percentage of GRHL2 enhancers with an overlapping CtIP enhancer: ", round(overlaps_grhl2_perc), "%"))
  
}



##


## Distances among CtIP and GRHL2 summits

# TODO: implement
# For those CtIP and GRHL2 enhancers overlapping (see above), save the HitsObject (instead of vector of overlapping regions)
# Store the summit of the overlapping enahncers (both CtIP and GRHL2)
# Plot a distribution of the distances among the summits of overlapping enhancers


# Alternative 
for(clust in 1:7){
  print(paste0("Checking overlaps in cluster: ", clust))
  
  # Analyze one cluster at a time to see how much the overlap across markers changes 
  enh_all_clust <- enh_all %>% lapply(., function(x){x %>% filter(cluster_new == clust)})
  left_join(enh_all_clust$CtIP, enh_all_clust$GRHL2[, c("chrom", "summit")], by = "chrom", 
            relationship = "many-to-many") # for each CtIP-enh on one chrom, we expect multiple GRHL2-enh on the same chrom 
  
  # TODO: implement
  # Store joined df 
  # Group joined df by "name" (chr:summit CtIP)
  # For each group (CtIP-enh), compute the distance of CtIP-enh summit from each GRHL2-enh summit found on the same chrom
  # For each group (CtIP-enh), pick the 'matching' GRHL2-enh summit CLOSEST to the CtIP-enh summit (min dist)
  # Plot the distribution of distances or some stats 
  # Given the selected 'matching' GRHL2 enhancers, check how many GRHL2 enhancers do not appear as close-enough (similar to non-overlapping)
  
  # Could it be usefult to do it also in the opposite direction? GRHL2 vs. CtIP?
  
}








