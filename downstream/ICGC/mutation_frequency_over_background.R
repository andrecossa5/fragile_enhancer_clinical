
library(tidyverse)
library(GenomicRanges)

SEED <- 4321
set.seed(SEED)

path_enhancers <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/")
path_SSMs <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.with_AFs.tsv")
path_chrom_sizes <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/hg19.chrom.txt")

MARKERS <- c("CtIP", "GRHL")
WIN <- 50
Mb <- 1000000


##


# Read enhancers 
columns_names <- c("chrom", "start", "end", "cluster")
enh_all <- list("CtIP" = read_tsv(fs::path(path_enhancers, "Cluster_CtIP_Enh_All.txt"), col_names = columns_names, comment = "#"), 
                "GRHL" = read_tsv(fs::path(path_enhancers, "Cluster_GRHL_Enh_All.txt"), col_names = columns_names, comment = "#"))
# add summit & name
enh_all <- lapply(enh_all, function(df){
  df$summit <- df$end
  df$name <- str_c(df$chrom, df$summit, sep=":")
  return(df)})


# Read variants
SSMs <- read_tsv(path_SSMs)
colnames(SSMs)[c(7,8)] <- c("start", "end")  
SSMs$chromosome <- paste("chr", SSMs$chromosome, sep = "")
SSMs <- SSMs[SSMs$end - SSMs$start == 0, ]

# Add mutatation ID equal to hartwig mutation ID
SSMs$ID <- str_sub(SSMs$chromosome, start = 4) %>% 
  str_c(., SSMs$end, sep = ":") %>% 
  str_c(., SSMs$reference_genome_allele, sep = "_") %>%
  str_c(., SSMs$mutated_to_allele, sep = "/")

SSMs_gr <- makeGRangesFromDataFrame(SSMs, keep.extra.columns = T)


# Read hg19 chrom sizes
standard_chrom <- paste0("chr", c(seq(1:22), "X", "Y"))
hg19_chrom_size <- read_tsv(path_chrom_sizes, col_names = c("chrom", "end")) %>% 
  dplyr::filter(., chrom %in% standard_chrom) %>%
  mutate(., start = 1) %>% relocate(., start, .before = end) 
hg19_chrom_size_gr <- makeGRangesFromDataFrame(hg19_chrom_size)


##

sample_var <- "icgc_sample_id"
fc_all_samples_markers <- list("CtIP", "GRHL")

for(marker in MARKERS){
  
  marker_enh <- enh_all[[marker]]
  marker_enh$start <- marker_enh$summit - WIN
  marker_enh$end <- marker_enh$summit + WIN
  enh_gr <- makeGRangesFromDataFrame(marker_enh, keep.extra.columns = T)
  
  # Store FoldChange(mutation rate over background)
  fc_all_samples <- c()
  
  for(sample in unique(SSMs[[sample_var]])){
    print(paste0("Analyzing sample: ", sample))
    
    # Get SNVs from 1 sample only 
    sample_SSMs_gr <- SSMs_gr[mcols(SSMs_gr)[[sample_var]] == sample]
    
    #hits_background
    hits_background <- findOverlaps(query=hg19_chrom_size_gr, subject=sample_SSMs_gr)
    back_SSMs <- cbind(data.frame(hg19_chrom_size_gr[queryHits(hits_background)]), data.frame(SSMs_gr[subjectHits(hits_background)]))
    tot_back_variants <- length(unique(back_SSMs$ID))
    tot_back_length_Mb <- round(sum(hg19_chrom_size$end) / Mb)
    # Compute background mutation frequency x Mb: number of SNVs x Mb (whole genome)     
    back_mut_rate <- round(tot_back_variants / tot_back_length_Mp, 3)
    
    #hits_enh
    hits_enh <- findOverlaps(query=enh_gr, subject=sample_SSMs_gr)
    print(paste0("Number of hits found for enhancers: ", length(hits_enh)))
    enh_SSMs <- cbind(data.frame(enh_gr[queryHits(hits_enh)]), data.frame(SSMs_gr[subjectHits(hits_enh)]))
    tot_enh_variants <- length(unique(enh_SSMs$ID))
    tot_enh_length_Mb <- round(sum(marker_enh$end - marker_enh$start) / Mb)
    # Compute enhancer mutation frequency x Mb: number of SNVs x Mb (enhancers)     
    enh_mut_rate <- round(tot_enh_variants / tot_enh_length_Mb, 3)
    
    # Fold change ( enh_mut_rate / back_mut_rate)
    fc <- enh_mut_rate / back_mut_rate
    fc_all_samples <- c(fc_all_samples, fc)
  }
  fc_all_samples_markers[[marker]] <- fc_all_samples
}


##


# Plot FC distribution for all samples
for(marker in MARKERS){
  df <- data.frame(fc = fc_all_samples_markers[[marker]]) 
  p <- df %>%
    ggplot(., aes(x = 1,y = fc))+
    geom_jitter(width = 0.3)+
    scale_x_continuous(limits = c(0,4), labels = c("", "enh", "", "", ""))+
    theme_light()+
    geom_hline(yintercept = 1, col = "grey", linetype = 2)+
    labs(x = "", y = "Fold-change over background", 
         title = paste0("Mutation frequency in ", marker, " enhancers"))
  print(p)
  
  # Same jitter_plot, but in log2 scale
  df_log <- data.frame(log2_fc = log2(fc_all_samples_markers[[marker]])+0.1) # adding a pseudocount to account for 0/X
  p_log <- df_log %>%
    ggplot(., aes(x = 1,y = log2_fc))+
    geom_jitter(width = 0.3)+
    scale_x_continuous(limits = c(0,4), labels = c("", "enh", "", "", ""))+
    scale_y_continuous(limits = c(-4, 4))+
    theme_light()+
    geom_hline(yintercept = 0, col = "grey", linetype = 2)
  print(p_log)
}





##


# TODO: 
# - Add distribution of mutation frequency for random regions
# - Plot distribution for high vs. low enhancers 
# - Plot distribution x cluster of enhancers 

# - What about WIN = 50 bp? - Check multiple windows 
# - Exclude hypermutated regions under physiological conditions?

# - implement for hartwig as well
# - implement for hartwig vs. ICGC

# - Define controls (what about inactive GRHL2/CtIP enhancers?)

