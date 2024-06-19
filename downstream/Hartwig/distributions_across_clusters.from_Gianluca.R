
library(tidyverse)
library(ggplot2)

path_overlaps <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/")
path_results_data <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/")

MARKERS <- c("CtIP", "GRHL")
save_dist <- F


##


# Extract distances of all PCAWG mutations from enhancer summits
g.icgc.overlaps <- list(
  "CtIP" = read_tsv(fs::path(path_overlaps, "Table_enh_SSMs_CtIP.all_overlaps.3kb_WIN.from_Gianluca.tsv")), 
  "GRHL" = read_tsv(fs::path(path_overlaps, "Table_enh_SSMs_GRHL.all_overlaps.3kb_WIN.from_Gianluca.tsv"))
)

for(marker in MARKERS){
  print(paste0("Computing distribution of SNVs distnaces from ", marker, " enhancers summit"))    
  all_distances <- data.frame(dist = as.numeric(g.icgc.overlaps[[marker]]$start_sbj) - as.numeric(g.icgc.overlaps[[marker]]$summit), 
                              cluster = g.icgc.overlaps[[marker]]$cluster)
  if(save_dist == T){
    all_distances %>% write_tsv(., fs::path(path_results_data, paste0("distances.snvs_distribution_over_enhancers.", marker, ".all_clusters.from_gianluca.tsv")))}
}


##


# Plot density of variants over enhancer regions
hart.ctip.dist <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/data/distances.snvs_distribution_over_enhancers.CtIP.all_clusters.tsv")
hart.grhl.dist <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/data/distances.snvs_distribution_over_enhancers.GRHL.all_clusters.tsv")
icgc.ctip.dist <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/distances.snvs_distribution_over_enhancers.CtIP.all_clusters.from_gianluca.tsv")
icgc.grhl.dist <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/distances.snvs_distribution_over_enhancers.GRHL.all_clusters.from_gianluca.tsv")
all_dfs <- list("CtIP" = list("hart" = hart.ctip.dist, "icgc" = icgc.ctip.dist), 
                "GRHL" = list("hart" = hart.grhl.dist, "icgc" = icgc.grhl.dist))

max_length <- max(sapply(list(hart.ctip.dist, hart.grhl.dist, icgc.ctip.dist, icgc.grhl.dist), function(x){dim(x)[1]}))

# Densities hartwig vs. icgc - ALL enhancers
df_full <- data.frame(
  "hart.ctip" = c(hart.ctip.dist$dist, rep(NA, max_length-dim(hart.ctip.dist)[1])), 
  "hart.grhl" = c(hart.grhl.dist$dist, rep(NA, max_length-dim(hart.grhl.dist)[1])),
  "icgc.ctip" = c(icgc.ctip.dist$dist, rep(NA, max_length-dim(icgc.ctip.dist)[1])),
  "icgc.grhl" = c(icgc.grhl.dist$dist, rep(NA, max_length-dim(icgc.grhl.dist)[1]))
)

for(marker in tolower(MARKERS)){
  df <- df_full[, endsWith(colnames(df_full), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("All clusters of enhancers - ", marker))
  print(d)
  
  #ggsave(plot=d, filename = paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/plots/density.snvs_distribution_over_enhancers.", marker,".all_clusters.hartwig_vs_icgc.png"), device = "png", width = 7, height = 5)
  
}


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

# x cluster high-low
df_high <- data.frame(
  "hart.ctip" = c(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$high, ]$dist, rep(NA, max_length-dim(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$high, ])[1])), 
  "hart.grhl" = c(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$high, ]$dist, rep(NA, max_length-dim(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$high, ])[1])),
  "icgc.ctip" = c(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$high, ]$dist, rep(NA, max_length-dim(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$high, ])[1])),
  "icgc.grhl" = c(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$high, ]$dist, rep(NA, max_length-dim(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$high, ])[1]))
)
df_low <- data.frame(
  "hart.ctip" = c(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$low, ]$dist, rep(NA, max_length-dim(hart.ctip.dist[hart.ctip.dist$cluster %in% cluster_groups$CtIP$low, ])[1])), 
  "hart.grhl" = c(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$low, ]$dist, rep(NA, max_length-dim(hart.grhl.dist[hart.grhl.dist$cluster %in% cluster_groups$GRHL$low, ])[1])),
  "icgc.ctip" = c(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$low, ]$dist, rep(NA, max_length-dim(icgc.ctip.dist[icgc.ctip.dist$cluster %in% cluster_groups$CtIP$low, ])[1])),
  "icgc.grhl" = c(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$low, ]$dist, rep(NA, max_length-dim(icgc.grhl.dist[icgc.grhl.dist$cluster %in% cluster_groups$GRHL$low, ])[1]))
)

# High 
for(marker in tolower(MARKERS)){
  df <- df_high[, endsWith(colnames(df_high), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("High clusters - ", marker))
  print(d)
  
  #ggsave(plot=d, filename = paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/plots/density.snvs_distribution_over_enhancers.", marker,".high_cluster.hartwig_vs_icgc.png"), device = "png", width = 7, height = 5)
}

# Low
for(marker in tolower(MARKERS)){
  df <- df_low[, endsWith(colnames(df_low), marker)]
  d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
    ggplot(., aes(x = distances, color = source, fill = source))+
    geom_density(bw = 100, alpha = 0.2)+
    theme_light()+
    labs(title = paste0("Low clusters - ", marker))
  print(d)
  
  #ggsave(plot=d, filename = paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/plots/density.snvs_distribution_over_enhancers.", marker,".high_cluster.hartwig_vs_icgc.png"), device = "png", width = 7, height = 5)
}


##


# Each cluster 

plot_list <- list(
  "CtIP" = setNames(vector(mode="list", length=length(unlist(cluster_groups[["CtIP"]]))), unlist(cluster_groups[["CtIP"]])),
  "GRHL" = setNames(vector(mode="list", length=length(unlist(cluster_groups[["GRHL"]]))), unlist(cluster_groups[["GRHL"]]))
)

for(marker in MARKERS){
  clusts <- unlist(cluster_groups[[marker]])
  for(clust in clusts){
    print(clust)
    df1 <- all_dfs[[marker]][["hart"]] %>% filter(cluster == clust)
    df2 <- all_dfs[[marker]][["icgc"]] %>% filter(cluster == clust)
    
    max_length <- max(dim(df1)[1], dim(df2)[1])
    df <- data.frame("hart" = c(df1$dist, rep(NA, max_length-length(df1$dist))), 
                     "icgc" = c(df2$dist, rep(NA, max_length-length(df2$dist))))
    
    d <- df %>% pivot_longer(everything(), names_to = "source", values_to = "distances") %>%
      ggplot(., aes(x=distances, group = source, color = source))+
      geom_density(bw=100, alpha=0.1)+
      theme_light()+
      labs(
        title = paste0(clust, " ", marker), 
        x = "distance from summit", y = "")
    print(d)
    
    plot_list[[marker]][[clust]] <- d
    }
}

for(marker in (MARKERS)){
  cowplot::plot_grid("title", plotlist = plot_list[[marker]], ncol = 1) %>%
    ggsave(., filename = paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/snvs_distribution_over_enhancers/plots/density.snvs_distribution_over_enhancers.", marker,".each_cluster.hartwig_vs_icgc.from_Gianluca.png"), device = "png", width = 7, height = 22)
}


