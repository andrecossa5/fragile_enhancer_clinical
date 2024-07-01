library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)

path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/")
  
SEED <- 4321
set.seed(SEED)


##

AF <- F

##


all.vcf.entries <- as.data.frame(matrix(nrow = 0, ncol = 13))
cols <- c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "TYPE", "PARID", "ID")
colnames(all.vcf.entries) <- cols

for(dir in list.dirs(path_input_main, recursive = F)){
  dir_name <- str_split(dir, pattern = "/", simplify = T)[length(str_split(dir, pattern = "/", simplify = T))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Read VCF file:
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.sv.vcf.gz")
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Keep only variants passing quality filters ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Extract somatic mutations ranges 
  somatic.vcf.filt.gr <- somatic.vcf.filt@rowRanges
  seqlevels(somatic.vcf.filt.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.gr))
  
  # Get column indicating breakpoints Partner ID
  if("PARID" %in% colnames(info(somatic.vcf.filt))){
    partner_id <- "PARID"
    partner_id_df <- data.frame("PARID" = info(somatic.vcf.filt)[, partner_id])
  } else {
    partner_id <- "MATEID"
    partner_id_df <- data.frame( "PARID" = sapply(info(somatic.vcf.filt)[, partner_id], function(x) if (length(x) == 0) NA else unlist(x)))
  }
  
  df_meta <- cbind(fixed(somatic.vcf.filt), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.gr)), 
                              "TYPE" = info(somatic.vcf.filt)[, 'EVENTTYPE']))
  
  df_meta <- cbind(df_meta, partner_id_df)
  df_meta$ID <- names(somatic.vcf.filt.gr)
  
  # For some reason, AF has 2 columns
  if(AF == T){
    df_AF <- lapply(info(somatic.vcf.filt)[, 'PURPLE_AF'], function(x) as.data.frame(t(x))) 
    df_AF <- do.call(rbind, df_AF)
    colnames(df_AF) <- c("PURPLE_AF_1", "PURPLE_AF_2")
    df_meta <- cbind(df_meta, df_AF)
  }
  mcols(somatic.vcf.filt.gr) <- df_meta
  
  somatic.vcf.filt.df <- data.frame(somatic.vcf.filt.gr)
  all.vcf.entries <- rbind(all.vcf.entries, somatic.vcf.filt.df)
}

# Save 
all.vcf.entries %>% write_tsv(., fs::path(path_results, paste0("Hartwig_all_stsms_info.tsv")))
