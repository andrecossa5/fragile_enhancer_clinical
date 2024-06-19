library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)

path_input_main <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/data")
path_results <- fs::path("/hpcnfs/scratch/P_PGP_FRAGILE_ENHANCERS/results/data/")
path_icgc_ssms <- fs::path("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/genomics/pre_processed_ICGC/simple_somatic_mutation.open.matching_calls.tsv")
  
SEED <- 4321
set.seed(SEED)
stats_icgc <- T

##


## HARTWIG 

all.enh.vcf.stats <- c()

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
  
  # Save n. of variants x sample
  n_var <- length(somatic.vcf.filt.snvs)
  all.enh.vcf.stats <- c(all.enh.vcf.stats, n_var) 
}

# Save overlaps x marker 
data.frame("n_vars" = all.enh.vcf.stats) %>% write_tsv(., fs::path(path_results, paste0("Stats_general.Hartwig_snvs_x_sample.tsv")))


##


## ICGC

# FIXME: filter only SNVs
if(stats_icgc ==T){
  icgc_ssms <- read_tsv(path_icgc_ssms)
  n_vars_icgc <- icgc_ssms %>% filter(nchar(mutated_to_allele) == 1) %>% 
    group_by(icgc_donor_id) %>%
    dplyr::summarise(n_vars = length(unique(icgc_mutation_id)))
  n_vars_icgc <- n_vars_icgc[, "n_vars"]
}



##

n_vars_hart <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/Hartwig/others/Stats_general.Hartwig_snvs_x_sample.tsv")

summary(n_vars_icgc)
summary(n_vars_hart)

n_vars_icgc$group <- "icgc"
n_vars_hart$group <- "hart"
rbind(n_vars_icgc, n_vars_hart) %>% 
  ggplot(., aes(x = n_vars, fill=group))+
  scale_x_log10()+
  geom_boxplot()+
  coord_flip()+
  theme_light()
