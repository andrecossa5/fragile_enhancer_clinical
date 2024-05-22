
#

library(tidyverse)
library(ggplot2)


path_overlaps <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data")

SEED <- 4321
set.seed(SEED) 
           
WIN <- 3000
MARKERS <- c("CtIP", "GRHL")


##


# Read enhancers x variants table overlaps 
overlaps <- setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
for(marker in MARKERS){
  overlaps[[marker]] <- read_tsv(fs::path(path_overlaps, paste0("Table_enh_SSMs_", marker, ".all_overlaps.", WIN/1000, "kb_WIN.tsv")))  
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

for(marker in MARKERS){
  # Distance of SNV from enhancers summit. start_sbj and end_abj are equal
  all_distances <- data.frame(dist = overlaps[[marker]]$start_sbj - overlaps[[marker]]$summit, 
                              cluster = overlaps[[marker]]$cluster)
  
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
  print(d_all)
  
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
    print(d_high_low)
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
    print(d_each)
  }
}



