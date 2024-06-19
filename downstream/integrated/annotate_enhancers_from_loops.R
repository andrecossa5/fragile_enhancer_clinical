
library(tidyverse)

SEED <- 4321
set.seed(SEED)


##

# All loops contain at least 1 enhancer #

for(kb in c(2,4)){
  print(paste0("Annotating enhancers from ", kb, " kb loops"))
  loops_table <-  read_tsv(sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb/data/tables/%skb_Unified_table.SCR_plus_KD_counts.all_anno_loops.ENH_DEGs_any.tsv",kb,kb)) %>%
    suppressMessages()
  
  # Select only loops of interest
  loi <- loops_table %>% 
    filter(is.na(counts.kd)) %>% #Select SCR-specific loops 
    filter((!is.na(name1) & !is.na(gene_name2) & DE2 == "Down") | (!is.na(name2) & !is.na(gene_name1) & DE1 == "Down")) #Select ENH-DOWN loops only
  
  enh_bin1 <- loi[!is.na(loi$gene_name2), ]$name1 
  enh_bin2 <- loi[!is.na(loi$gene_name1), ]$name2
  enh_oi <- c(enh_bin1[!is.na(enh_bin1)], enh_bin2[!is.na(enh_bin2)])
  enh_oi <- data.frame("name" = unique(enh_oi))
  
  print(paste0("Total number of enhancers from SCR-specific-Down loops, ", kb, " kb: ", dim(enh_oi)[1]))
  
  path_output <- sprintf("/Users/ieo6983/Desktop/enhancers_project/Analyses/loops/results/%skb/data/anno_enhancers/", kb)
  file_name <- sprintf("%skb_GRHL2_enhancers.from_SCR_specific_loops.linked_to_DOWN_DEGs.tsv", kb)
  #enh_oi %>% write_tsv(., paste0(path_output, file_name))
}


