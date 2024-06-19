
### Compute stats ###

library(tidyverse)

path_overlaps_hart <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data")
path_overlaps_icgc <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/results/ICGC/enhancers_SSMs_overlaps/data/")

path_enhancers_ctip <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/others/Cluster_CtIP_Enh_All.txt")
path_enhancers_grhl <- fs::path("/hpcnfs/scratch/PGP/Ciacci_et_al/data/functional_genomics/others/Cluster_GRHL_Enh_All.txt")

path_results_data <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data")
path_results_plots <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/plots")

SEED <- 4321
set.seed(SEED) 

MARKERS <- c("CtIP", "GRHL")
SOURCES <- c("ICGC", "HART")
N_DON <- setNames(list(793, 1059), SOURCES) # hart = n. of folders 
N_ENH <- list("CtIP" = dim(read_tsv(path_enhancers_ctip,col_names=F,comment="#"))[1], 
              "GRHL" = dim(read_tsv(path_enhancers_grhl,col_names=F,comment="#"))[1])

save_plots <- F
# TODO: add code to save tables with stats when needed
save_tables <- F 


##


# Read enhancers x variants table overlaps 
overlaps <- list(
  "ICGC" = setNames(vector(mode="list", length=length(MARKERS)), MARKERS), 
  "HART" = setNames(vector(mode="list", length=length(MARKERS)), MARKERS)
)

for(marker in MARKERS){
  column_names <- c(
    paste0(c("seqnames", "start", "end", "width", "strand", "cluster", "summit", "name"), "_enh"), 
    paste0(c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "SAMPLE", "ID"), "_snv")
  )
  
  # ICGC
  overlaps$ICGC[[marker]] <- read_tsv(fs::path(path_overlaps_icgc, paste0("Table_enh_SSMs_", marker, ".all_overlaps.3kb_WIN.tsv"))) %>% 
    dplyr::select(c(1:13,21:22, 15, 14)) 
  colnames(overlaps$ICGC[[marker]]) <- column_names
  
  # HART
  overlaps$HART[[marker]] <- read_tsv(fs::path(path_overlaps_icgc, paste0("Table_enh_SSMs_", marker, ".all_overlaps.3kb_WIN.tsv"))) %>% 
    dplyr::select(c(1:13, 14:15, 18:19)) 
  colnames(overlaps$HART[[marker]]) <- column_names
}

# Add distances of SNVs from enh-summit
for(i in 1:length(overlaps)){
  for(marker in MARKERS){
    overlaps[[i]][[marker]]$dist <- abs(overlaps[[i]][[marker]]$start_snv - overlaps[[i]][[marker]]$summit_enh)
  }
}


##


df_perc <- data.frame(
  "win" = NA, 
  "marker" = NA,
  "source" = NA, 
  "perc_don" = NA, 
  "perc_enh" = NA
)

df_box1 <- list(
  "CtIP" = matrix(nrow=max(unlist(N_DON)),ncol=1), 
  "GRHL" = matrix(nrow=max(unlist(N_DON)),ncol=1)
)
col_names <- list("CtIP" = c(), "GRHL" = c())

df_box2 <- list(
  "CtIP" = matrix(nrow=max(unlist(N_ENH$CtIP)),ncol=1), 
  "GRHL" = matrix(nrow=max(unlist(N_ENH$GRHL)),ncol=1)
)
col_names2 <- list("CtIP" = c(), "GRHL" = c())

for(WIN in c(3000, 2000, 1000, 500, 100, 50)){
  cat("\n")
  print(paste0("Computing stats for window: ", WIN))
  
  for(marker in MARKERS){
    cat("\n")
    print(paste0("Marker = ", marker))
    
    for(source in SOURCES){
      df_perc_add <- c(WIN, marker, source)
      overlaps_filt <- overlaps[[source]][[marker]] %>% dplyr::filter(dist <= WIN)
      
      # Percentage of donors with at least 1 mutated enhancer
      n_don <- length(unique(overlaps_filt$SAMPLE_snv))
      perc <- round((n_don / N_DON[[source]]) * 100)
      print(paste0("Number of ", source, " donors with at least 1 ", marker, " enhancer mutated: ", n_don, " over ", N_DON[[source]], " (", perc, "%)"))
      df_perc_add <- c(df_perc_add, perc)
      
      # Percentage of enhancers mutated in at least 1 donor
      n_enh <- length(unique(overlaps_filt$name_enh))
      perc <- round((n_enh / N_ENH[[marker]]) * 100)
      print(paste0("Number of ", marker, " enhancers mutated in at least 1 ", source, " donor: ", n_enh, " over ", N_ENH[[marker]], " (", perc, "%)"))
      df_perc_add <- c(df_perc_add, perc)
      
      df_perc <- rbind(df_perc, df_perc_add)
      df_perc_add <- c()
      
      # Number of enhancers mutated x donor
      n_enh_dist <- overlaps_filt %>% group_by(SAMPLE_snv) %>%
        dplyr::summarise("n_enh_mut" = length(unique(name_enh)))
      n_enh_dist <- n_enh_dist$n_enh_mut
      # Add donors with 0 enhancers mutated
      n_enh_dist <- c(n_enh_dist, rep(NA, (N_DON[[source]] - length(n_enh_dist))))
      # Add distribution to list of matrices
      col_name <- paste(WIN, marker, source, "n_enh_mut", sep=".")
      col_names[[marker]] <- c(col_names[[marker]], col_name)
      df_box1[[marker]] <- cbind(df_box1[[marker]], n_enh_dist)
      
      # Number of mutations x enhancer x donors 
      n_muts_enh <- overlaps_filt %>% group_by(name_enh, SAMPLE_snv) %>% 
        dplyr::summarise(n_muts_enh_x_don = n())
      n_muts_enh_vec <- n_muts_enh$n_muts_enh_x_don  
      # Add enhancers with 0 enhancers mutated
      n_muts_enh_x_don_dist <- c(n_muts_enh_vec, rep(NA, (N_ENH[[marker]] - length(unique(n_muts_enh$name_enh)))))
      # Add distribution to list of matrices
      col_name2 <- paste(WIN, marker, source, "n_muts_enh_x_don", sep=".")
      col_names2[[marker]] <- c(col_names2[[marker]], col_name2)
      df_box2[[marker]] <- cbind(df_box2[[marker]], n_muts_enh_x_don_dist)
    }
  }
}


##


# Plot percentages 
# TODO: check if correct, since CtIP and GRHL appear exactly the same  
for(marker in MARKERS){
  df_perc <- na.omit(df_perc)
  df_perc$win <- factor(df_perc$win, levels = c(3000, 2000, 1000, 500, 100, 50))
  df_perc$source <- factor(df_perc$source, levels = c("ICGC", "HART"))
  df_perc$perc_don <- as.numeric(df_perc$perc_don)
  df_perc$perc_enh <- as.numeric(df_perc$perc_enh)
  
  # Plot percentage of donors with at least 1 mutated enhancer across datasets and windows
  p_don <- df_perc %>% dplyr::filter(marker == marker) %>% 
    ggplot(., aes(x=win, y=perc_don, group=source, fill=source))+
    geom_bar(stat = "identity", position = "dodge")+
    labs(title = "Percentage of donors with at least 1 mutated enhnacer", 
         subtitle = marker, 
         x = "window", y = "%")+
    theme_light()+
    scale_fill_manual(values = c("ICGC" = "#d8a6a6", "HART" = "red4"))
  #print(p_don)
  
  # Plot percentage of enhancers mutated in at least 1 donor across datasets and windows
  p_enh <- df_perc %>% dplyr::filter(marker == marker) %>% 
    ggplot(., aes(x=win, y=perc_enh, group=source, fill=source))+
    geom_bar(stat = "identity", position = "dodge")+
    labs(title = "Percentage of enhancers mutated in at least 1 donor", 
         subtitle = marker, 
         x = "window", y = "%")+
    theme_light()+
    scale_fill_manual(values = c("ICGC" = "#d8a6a6", "HART" = "red4"))
  #print(p_enh)
  
  if(save_plots == T){
    ggsave(plot=p_don, filename=fs::path(path_results_plots, paste0("bars.perc_donors_mut.", marker, ".all_clusters.all_windows.png")), 
           device = "png", width = 7.791667, height = 6.666667)
    ggsave(plot=p_enh, filename=fs::path(path_results_plots, paste0("bars.perc_enhancers_mut.", marker, ".all_clusters.all_windows.png")), 
           device = "png", width = 7.791667, height = 6.666667)
  }
}


# Plot Distributions
for(marker in MARKERS){
  df <- data.frame(df_box1[[marker]]) %>% select(-V1)
  colnames(df) <- col_names[[marker]]
  df_box1[[marker]] <- df
  rm(df)
  
  df <- data.frame(df_box2[[marker]]) %>% select(-V1)
  colnames(df) <- col_names2[[marker]]
  df_box2[[marker]] <- df
  rm(df)
}

for(marker in MARKERS){
  # Number of enhancers mutated across donors
  df1 <- data.frame("groups" = colnames(df_box1[[marker]]))
  df2 <- data.frame(str_split(string = colnames(df_box1$CtIP), pattern = "\\.", simplify = T))
  colnames(df2) <- c("win", "marker", "source", "name")
  df_legend <- cbind(df1, df2)
  
  df_long <- df_box1[[marker]] %>% pivot_longer(everything(), names_to = "groups", values_to = "value") %>%
    left_join(., df_legend, by = "groups")
  df_long$win <- factor(df_long$win, levels = c(3000, 2000, 1000, 500, 100, 50))
  df_long$source <- factor(df_long$source, levels = c("ICGC", "HART"))
    
  p1 <- df_long %>% ggplot(., aes(x = win, y = value, fill = source))+
    geom_boxplot()+
    labs(title = "Number of enhancers mutated x donor donors", 
         subtitle = marker, 
         x = "window", y = "")+
    theme_light()+
    scale_fill_manual(values = c("ICGC" = "#d8a6a6", "HART" = "red4"))
  #print(p1)
  
  # Number of mutations x enhancer x donor
  df1 <- data.frame("groups" = colnames(df_box2[[marker]]))
  df2 <- data.frame(str_split(string = colnames(df_box2$CtIP), pattern = "\\.", simplify = T))
  colnames(df2) <- c("win", "marker", "source", "name")
  df_legend <- cbind(df1, df2)
  
  df_long <- df_box2[[marker]] %>% pivot_longer(everything(), names_to = "groups", values_to = "value") %>%
    left_join(., df_legend, by = "groups")
  df_long$win <- factor(df_long$win, levels = c(3000, 2000, 1000, 500, 100, 50))
  df_long$source <- factor(df_long$source, levels = c("ICGC", "HART"))
  
  
  fac <- 1000000
  pseudoc <- 1
  p2 <- df_long %>% ggplot(., aes(x = win, y = value, fill = source))+
    geom_violin()+
    labs(title = "Number of mutations x enhancer x donor", 
         subtitle = marker, 
         x = "window", y = "")+
    theme_light()
  #+scale_fill_manual(c("ICGC" = "#d8a6a6", "HART" = "red4"))
  #print(p2)
 
  if(save_plots == T){
    ggsave(plot=p1, filename=fs::path(path_results_plots, paste0("box.n_enhancers_mut_x_don.", marker, ".all_clusters.all_windows.png")), 
           device = "png", width = 8.666667, height = 10.347222)
    ggsave(plot=p2, filename=fs::path(path_results_plots, paste0("violin.n_muts_x_enhancer_x_don.", marker, ".all_clusters.all_windows.png")), 
           device = "png", width = 8.666667, height = 10.347222)
  } 
}


##



