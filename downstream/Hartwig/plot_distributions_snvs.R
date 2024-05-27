
#

library(tidyverse)
library(ggplot2)


path_overlaps <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data")
path_results_data <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data")
path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/plots")
  
SEED <- 4321
set.seed(SEED) 
           
WIN <- 3000
MARKERS <- c("CtIP", "GRHL")


##


# Read enhancers x variants table overlaps 
overlaps <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)

# FIXME: 
# - overlaps file has duplicated col_names (seqnames, start, end, etc.). 
# - change find_overlaps_snvs.R to save the output with non-duplicated col_names
# - change code below: remove addition of column_names 

for(marker in MARKERS){
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER", "SAMPLE", "ID"), "_snv")
  )
  overlaps[[marker]] <- read_tsv(fs::path(path_overlaps, paste0(marker, "_enh.hartwig_snvs.overlap.WIN_", WIN, ".tsv")), col_names = column_names)  
}

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

                     
##


# Plot distirbution of SNVs over enhancers regions

high_low <- T
each_cluster <- T
save_plots <- F
save_dist <- T
  
for(marker in MARKERS){
  # Distance of SNV from enhancers summit. start_sbj and end_abj are equal
  cat("\n")
  print(paste0("Computing distribution of SNVs distnaces from ", marker, " enhancers summit"))
  all_distances <- data.frame(dist = as.numeric(overlaps[[marker]]$start_snv) - as.numeric(overlaps[[marker]]$summit_enh), 
                              cluster = overlaps[[marker]]$cluster_enh)
  if(save_dist == T){
    all_distances %>% write_tsv(., fs::path(path_results_data, paste0("distances.snvs_distribution_over_enhancers.", marker, ".all_clusters.tsv")))}
  
  
  # Plot density of SNVs distribution over enhancer region - ALL enhancers
  bw <- 100   
  d_all <- all_distances %>% ggplot(., aes(dist))+
    geom_density(fill = "lightblue3", color = "lightblue3",alpha = 0.7, bw = bw)+
    theme_light()+
    scale_x_continuous(breaks = seq(-3000, 3000, 500))+
    labs(
      title = paste0("SNVs distribution over ", marker, " enhancers"),
      x = "distance from summit", y = ""
    )
  print("Plotting SNVs distribution over enhancer region - ALL enhancers")
  #print(d_all)
  
  if(save_plots == T){
  ggsave(plot=d_all, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".all_clusters.png")), 
         device = "png", width = 9.805556, height = 7.027778)}
  

  # Plot density of SNVs distribution over enhancer region - HIGH vs. LOW enhancers
  if(high_low ==T){
    dist_high <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$high) 
    dist_low <- all_distances %>% filter(cluster %in% cluster_groups[[marker]]$low) 
    
    max_dist_length <- max(dim(dist_high)[1], dim(dist_low)[1])
    df_grouped <- data.frame("high" = c(dist_high$dist, rep(NA, max_dist_length - dim(dist_high)[1])), 
                             "low" = c(dist_low$dist, rep(NA, max_dist_length - dim(dist_low)[1])))
    
    d_high_low <- df_grouped %>% pivot_longer(everything(), names_to = "cluster", values_to = "distances") %>%
      ggplot(., aes(x = distances, color = cluster, fill = cluster))+
      geom_density(bw = bw, alpha = 0.2)+
      theme_light()+
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
      )
    print("Plotting SNVs distribution over enhancer region - HIGH vs. LOW enhancers")
    #print(d_high_low)
    
    if(save_plots == T){
    ggsave(plot=d_high_low, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".high_low.png")), 
           device = "png", width = 9.805556, height = 7.027778)}
  }
  
  # Plot density of SNVs distribution over enhancer region - EACH cluster of enhancers
  if(each_cluster == T){
    d_each <- all_distances %>% 
      ggplot(., aes(x = dist, gorup = cluster, color = cluster))+
      geom_density(bw = bw, alpha =0.1)+
      theme_light()+
      labs(
        title = paste0("SNVs distribution over ", marker, " enhancers"), 
        x = "distance from summit", y = ""
        )
    print("Plotting SNVs distribution over enhancer region - EACH cluster of enhancers")
    #print(d_each)
    
    if(save_plots == T){
    ggsave(plot=d_each, filename=fs::path(path_results_plots, paste0("density.snvs_distribution_over_enhancers.", marker, ".each_cluster.png")), 
           device = "png", width = 9.805556, height = 7.027778)}
  }
}


# TODO: convert into a function 
# TODO: chnage script to accept arguments from command line

