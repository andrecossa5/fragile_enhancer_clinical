
library(fs)  # File manipulations
library(tidyverse)
library(GenomicRanges)
library(viridis)
library(gridExtra)

source("./Desktop/enhancers_project/Analyses/loops/loops_functions.R")

### Hi-ChIP Loops ###
kb <- 2

path_main <- "/Users/ieo6983/Desktop/fragile_enhancer_clinical"
path_data <- fs::path(path_main, "data") 
path_results <- fs::path(path_main, "results")
path_hichip <- fs::path(path_data, "functional_genomics/HiChip/filtered_loops/")
path_enhancers <- fs::path(path_data, "functional_genomics/others")
path_degs <- "/Users/ieo6983/Desktop/expression/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv"

path_main_input <- sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb/data", kb)
path_output <- fs::path(sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb", kb))


##


# Read loops
scr_enh_degs_only <- read_tsv(fs::path(path_main_input, sprintf("%skb_SCR.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv", kb)))
kd_enh_degs_only <- read_tsv(fs::path(path_main_input, sprintf("%skb_KD.anno_loops.GRHL2_enh_DEGs_prom_ONLY.tsv", kb)))

# condition-specific
scr_specific <- read_tsv(fs::path(path_main_input, sprintf("%skb_SCR_specific_loops.tsv", kb)))
kd_specific <- read_tsv(fs::path(path_main_input, sprintf("%skb_KD_specific_loops.tsv", kb)))

# Read GRHL2-bound enhancers
columns_names <- c("chrom", "start", "end", "cluster")
enh_grhl2 <- read_tsv(fs::path(path_enhancers, "Cluster_GRHL_Enh_All.txt"), col_names = columns_names, comment = "#")
# Pre-process: add summit & enhancer name 
enh_grhl2$summit <- enh_grhl2$end
enh_grhl2$name <- str_c(enh_grhl2$chrom, enh_grhl2$summit, sep=":")

# Add high-low category
clusters_high <- c("GRHL_cluster_1_Enh", "GRHL_cluster_2_Enh", "GRHL_cluster_4_Enh")
clusters_low <- c("GRHL_cluster_5.0", "GRHL_cluster_5.1_Enh", "GRHL_cluster_5.2_Enh", "GRHL_cluster_5.3_Enh")
enh_grhl2$group <- "high"
enh_grhl2[enh_grhl2$cluster %in% clusters_low, ]$group <- "low"


##


# Unbundle ambiguous overlaps
scr_enh_degs_only_filt <- unbundle_ambiguous_bins(scr_enh_degs_only)
kd_enh_degs_only_filt <- unbundle_ambiguous_bins(kd_enh_degs_only)


##


enh_degs_only <-  list("scr" = scr_enh_degs_only_filt, "kd" = kd_enh_degs_only_filt)

# SCR-specific enhancers connected to DEGs - N. assocaited to Up/Down regulated genes
t1 <- table(c(scr_enh_degs_only_filt$DE1, scr_enh_degs_only_filt$DE2)[!is.na(c(scr_enh_degs_only_filt$DE1, scr_enh_degs_only_filt$DE2))])

# KD-specific enhancers connected to DEGs
t2 <- table(c(kd_enh_degs_only_filt$DE1, kd_enh_degs_only_filt$DE2)[!is.na(c(kd_enh_degs_only_filt$DE1, kd_enh_degs_only_filt$DE2))])

loop_degs <- data.frame(group = c("scr", "kd"),
                        down = c(t1["Down"], t2["Down"]),
                        up = c(t1["Up"], t2["Up"]))
loop_degs$group <- factor(loop_degs$group, levels = c("scr", "kd"))

loop_degs %>% pivot_longer(-group) %>%
  ggplot(., aes(x=group, y = value, fill = name))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Loops association to DEGs", 
       x = "", y = "number of loops")+
  scale_fill_manual(values = c("purple3", "orange3"))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  scale_y_continuous(n.breaks = 10)


# DEGs involved in loops
DEGs <- read_tsv(path_degs)

degs_in_loops <- list()

for(cond in c("scr", "kd")){
  df_loops <- enh_degs_only[[cond]]

  degs_uniq <- data.frame("gene_name" = unique(c(na.omit(df_loops$gene_name1), na.omit(df_loops$gene_name2))))
  degs_uniq <- left_join(degs_uniq, DEGs[, c("gene_name", "DE", "log2FoldChange", "padj")], by = "gene_name")
  
  degs_in_loops[[cond]] <- degs_uniq
  #degs_uniq %>% write_tsv(., fs::path(path_output, sprintf("/data/%skb_%s.DEGs_in_enh_deg_loops",kb, toupper(cond)), ext = "tsv"))
}

# Number of DEGs involved in loops
for(cond in c("scr", "kd")){
  dfs_cond <- enh_degs_only[[cond]]
  
  degs_n_uniq <- length(unique(c(na.omit(dfs_cond$gene_name1), na.omit(dfs_cond$gene_name2))))
  print(sprintf("%s: %s", cond, degs_n_uniq))
  
  up1 <- dfs_cond %>% filter(DE1 == "Up")
  up2 <- dfs_cond %>% filter(DE2 == "Up")
  degs_up <- length(unique(na.omit(c(up1$gene_name1, up2$gene_name2))))
  down1 <- dfs_cond %>% filter(DE1 == "Down")
  down2 <- dfs_cond %>% filter(DE2 == "Down")
  degs_down <- length(unique(na.omit(c(down1$gene_name1, down2$gene_name2))))
  print(sprintf("%s - Down: %s; Up: %s", cond, degs_down, degs_up))
}


# Cluster composition 
for(cond in c("scr", "kd")){
  df_loops <- enh_degs_only[[cond]]
  clusters_in_loops <- c(df_loops$cluster1, df_loops$cluster2)[!is.na(c(df_loops$cluster1, df_loops$cluster2))]
  
  print(sprintf("Cluster composition of ENH-DEG loops - %s", toupper(cond)))
  print(table(clusters_in_loops))
  cat("\n")
  
  h <- as.data.frame(table(clusters_in_loops)) %>% 
    filter(clusters_in_loops %in% clusters_high) %>% 
    select(Freq) %>% colSums()
  l <- as.data.frame(table(clusters_in_loops)) %>% 
    filter(clusters_in_loops %in% clusters_low) %>% 
    select(Freq) %>% colSums()
  
  hr <- h / sum(enh_grhl2$group == "high")
  lr <- l / sum(enh_grhl2$group == "low")
  
  print(sprintf("Number of enh. in clusters HIGH: %s (%s) vs. number of enh. in clusters LOW: %s (%s)", 
                h, round(hr,3), 
                l, round(lr,3)))
  cat("\n")
  
  # Define observed counts
  observed_counts <- matrix(c(h, l, sum(enh_grhl2$group == "high")-h, sum(enh_grhl2$group == "low")-l), nrow = 2, byrow = TRUE)
  colnames(observed_counts) <- c("High", "Low")
  rownames(observed_counts) <- c("In loops", "Not in loops")
  chi_squared_test <- chisq.test(observed_counts)
  
  print(sprintf("Chi-square test p-value: %s", round(chi_squared_test$p.value,3)))
  cat("\n")
}


##


# SCR-specific enhancers connected to DEGs
scr_enh_degs_only_spec <- inner_join(scr_enh_degs_only_filt, scr_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
t1 <- table(c(scr_enh_degs_only_spec$DE1, scr_enh_degs_only_spec$DE2)[!is.na(c(scr_enh_degs_only_spec$DE1, scr_enh_degs_only_spec$DE2))])

# KD-specific enhancers connected to DEGs
kd_enh_degs_only_spec <- inner_join(kd_enh_degs_only_filt, kd_specific, by = c('seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'))
t2 <- table(c(kd_enh_degs_only_spec$DE1, kd_enh_degs_only_spec$DE2)[!is.na(c(kd_enh_degs_only_spec$DE1, kd_enh_degs_only_spec$DE2))])

cond_spec_degs <- data.frame(group = c("scr", "kd"),
                             down = c(t1["Down"], t2["Down"]), 
                             up = c(t1["Up"], t2["Up"]))
cond_spec_degs$group <- factor(cond_spec_degs$group, levels = c("scr", "kd"))

cond_spec_degs %>% pivot_longer(-group) %>%
  ggplot(., aes(x=group, y = value, fill = name))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_light()+
  labs(title = "Condition-specific loops & their association to DEGs", 
       x = "", y = "number of loops")+
  scale_fill_manual(values = c("purple3", "orange3"))+
  theme(axis.text.x = element_text(size = 12, color = "black"))


# Number of DEGs involved in condition-sepcific loops
cond_spec <-  list("scr" = scr_enh_degs_only_spec, "kd" = kd_enh_degs_only_spec)

for(cond in c("scr", "kd")){
  dfs_cond <- cond_spec[[cond]]
  
  degs_n_uniq <- length(unique(c(na.omit(dfs_cond$gene_name1), na.omit(dfs_cond$gene_name2))))
  print(sprintf("%s: %s", cond, degs_n_uniq))
  
  up1 <- dfs_cond %>% filter(DE1 == "Up")
  up2 <- dfs_cond %>% filter(DE2 == "Up")
  degs_up <- length(unique(na.omit(c(up1$gene_name1, up2$gene_name2))))
  down1 <- dfs_cond %>% filter(DE1 == "Down")
  down2 <- dfs_cond %>% filter(DE2 == "Down")
  degs_down <- length(unique(na.omit(c(down1$gene_name1, down2$gene_name2))))
  print(sprintf("%s - Down: %s; Up: %s", cond, degs_down, degs_up))
}

# DEGs involved in loops
DEGs <- read_tsv("./Desktop/expression/DEGs/Df_DEGs.df_LFC_sig.padj_0.05.log2FC_1.Up_and_Down.tsv")

degs_in_loops_spec <- list()

for(cond in c("scr", "kd")){
  df_loops <- cond_spec[[cond]]
  
  degs_uniq <- data.frame("gene_name" = unique(c(na.omit(df_loops$gene_name1), na.omit(df_loops$gene_name2))))
  degs_uniq <- left_join(degs_uniq, DEGs[, c("gene_name", "DE", "log2FoldChange", "padj")], by = "gene_name")
  
  degs_in_loops_spec[[cond]] <- degs_uniq
  #degs_uniq %>% write_tsv(., fs::path(path_output, sprintf("/data/%skb_%s.DEGs_in_enh_deg_loops.cond_spec",kb, toupper(cond)), ext = "tsv"))
}


# Cluster composition
for(cond in c("scr", "kd")){
  df_loops <- cond_spec[[cond]]
  clusters_in_loops <- c(df_loops$cluster1, df_loops$cluster2)[!is.na(c(df_loops$cluster1, df_loops$cluster2))]
  
  print(sprintf("Cluster composition of condition-specific ENH-DEG loops - %s", toupper(cond)))
  print(table(clusters_in_loops))
  cat("\n")
  
  h <- as.data.frame(table(clusters_in_loops)) %>% 
    filter(clusters_in_loops %in% clusters_high) %>% 
    select(Freq) %>% colSums()
  l <- as.data.frame(table(clusters_in_loops)) %>% 
    filter(clusters_in_loops %in% clusters_low) %>% 
    select(Freq) %>% colSums()
  
  hr <- h / sum(enh_grhl2$group == "high")
  lr <- l / sum(enh_grhl2$group == "low")
  
  print(sprintf("Number of enh. in clusters HIGH: %s (%s) vs. number of enh. in clusters LOW: %s (%s)", 
                h, round(hr,3), 
                l, round(lr,3)))
  cat("\n")
  
  # Define observed counts
  observed_counts <- matrix(c(h, l, sum(enh_grhl2$group == "high")-h, sum(enh_grhl2$group == "low")-l), nrow = 2, byrow = TRUE)
  colnames(observed_counts) <- c("High", "Low")
  rownames(observed_counts) <- c("In loops", "Not in loops")
  chi_squared_test <- chisq.test(observed_counts)
  
  print(sprintf("Chi-square test p-value: %s", round(chi_squared_test$p.value,3)))
  cat("\n")
}


##


# GSEA on DEGs-associated loops
library(fgsea)

# MSigDB Gene Set
hallmark <- "/Users/ieo6983/Desktop/expression/GSEA/DB/hallmark_gene_sets.h.all.v2023.2.Hs.symbols.gmt"

## DEGs in ENH_DEG loops
gsea_degs <- list()

for(cond in c("scr", "kd")){
  # List of DEGs, ranked by log2FC
  ranked_degs_df <- degs_in_loops[[cond]] %>% arrange(., -log2FoldChange)
  ranked_degs_list <- ranked_degs_df$log2FoldChange
  names(ranked_degs_list) <- ranked_degs_df$gene_name
  
  # P-value threhsold to filter pathways
  p_val <- 0.05
  
  gsea_degs[[cond]] <- GSEA(gene_list = ranked_degs_list, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
}

gsea_degs$scr
gsea_degs$kd

# All DEGs (SCR & KD) in ENH-DEG loops
ranked_all_df <- rbind(degs_in_loops$scr, degs_in_loops$kd) %>% filter(!duplicated(.)) %>% arrange(., -log2FoldChange)
ranked_all <- ranked_all_df$log2FoldChange
names(ranked_all) <- ranked_all_df$gene_name
gsea_degs_all <- GSEA(gene_list = ranked_all, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
gsea_degs_all


## DEGs in ENH-DEG condition-specific loops
gsea_degs_spec <- list()

for(cond in c("scr", "kd")){
  # List of DEGs, ranked by log2FC
  ranked_degs_df <- degs_in_loops_spec[[cond]] %>% arrange(., -log2FoldChange)
  ranked_degs_list <- ranked_degs_df$log2FoldChange
  names(ranked_degs_list) <- ranked_degs_df$gene_name
  
  # P-value threhsold to filter pathways
  p_val <- 0.05
  
  gsea_degs_spec[[cond]] <- GSEA(gene_list = ranked_degs_list, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
}

gsea_degs_spec$scr
gsea_degs_spec$kd

# All DEGs (SCR & KD) in ENH-DEG condition-specific loops
ranked_all_df <- rbind(degs_in_loops$scr, degs_in_loops_spec$kd) %>% filter(!duplicated(.)) %>% arrange(., -log2FoldChange)
ranked_all <- ranked_all_df$log2FoldChange
names(ranked_all) <- ranked_all_df$gene_name
gsea_degs_all <- GSEA(gene_list = ranked_all, GO_file = hallmark, p_val = p_val, minSize = 0, maxSize = Inf, collapse = F)
gsea_degs_all


##



