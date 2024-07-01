library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)

path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/")
  
SEED <- 4321
set.seed(SEED)


##


all.vcf.entries <- as.data.frame(matrix(nrow = 0, ncol = 13))
cols <- c("seqnames", "start", "end", "width", "strand", "REF", "ALT", "QUAL", "FILTER",  "SAMPLE", "PURPLE_AF", "AF", "ID")
colnames(all.vcf.entries) <- cols

for(dir in list.dirs(path_input_main, recursive = F)){
  dir_name <- str_split(dir, pattern = "/", simplify = T)[length(str_split(dir, pattern = "/", simplify = T))]
  cat("\n")
  print(paste0("Iterating over sample: ", dir_name))
  
  # Read VCF file:
  somatic.filename <- paste0(dir, "/purple/", dir_name, ".purple.somatic.vcf.gz")
  somatic.vcf <- readVcf(somatic.filename, genome = "hg19")
  
  # Keep only variants passing quality filters ("PASS")
  qc_filt <- fixed(somatic.vcf)[, "FILTER"] == "PASS"
  somatic.vcf.filt <- somatic.vcf[qc_filt]
  print("-- Dropping variants NOT passing QC filters --")
  print(paste0("Initial number of variants: ", length(somatic.vcf), " | Number of variants after filtering: ", length(somatic.vcf.filt)))
  
  # Keep only SNVs
  snv_filt <- (nchar(fixed(somatic.vcf.filt)[, "REF"]) == 1) & unlist((nchar(fixed(somatic.vcf.filt)[, "ALT"]) == 1))
  somatic.vcf.filt.snvs <- somatic.vcf.filt[snv_filt]
  print("-- Retaining only SNVs --")
  print(paste0("Initial number of variants: ", length(somatic.vcf.filt), " | Number of variants after filtering: ", length(somatic.vcf.filt.snvs)))
  
  # Extract somatic mutations ranges 
  somatic.vcf.filt.snvs.gr <- somatic.vcf.filt.snvs@rowRanges
  seqlevels(somatic.vcf.filt.snvs.gr) <- paste0("chr", seqlevels(somatic.vcf.filt.snvs.gr))
  
  df_meta <- cbind(fixed(somatic.vcf.filt.snvs), 
                   data.frame("SAMPLE" = rep(dir_name, length(somatic.vcf.filt.snvs.gr))), 
                   data.frame("PURPLE_AF" = info(somatic.vcf.filt.snvs)[, 'PURPLE_AF']), 
                   data.frame("AF" = geno(somatic.vcf.filt.snvs)$AF[, 2])
  )
  df_meta$ID <- names(somatic.vcf.filt.snvs.gr)
  mcols(somatic.vcf.filt.snvs.gr) <- df_meta
  
  somatic.vcf.filt.snvs.df <- data.frame(somatic.vcf.filt.snvs.gr)
  all.vcf.entries <- rbind(all.vcf.entries, somatic.vcf.filt.snvs.df)
}

# Save 
all.vcf.entries %>% write_tsv(., fs::path(path_results, paste0("Hartwig_all_snvs_info.tsv")))
