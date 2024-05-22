
library(tidyverse)
library(rtracklayer)

### Gene annotation - USCS Table browser - ENSEMBL Genes ###

anno <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/Gene_anno_hg19_ENSEMBL.tsv")

# Transcripts x gene
tx_x_g <- anno %>% group_by(name2) %>% 
  dplyr::summarise(., n = n())
summary(tx_x_g)
hist(table(anno$name2))

# Get TSSs
TSS_plus <- anno[anno$strand == "+", c("name2", "chrom", "strand", "txStart")]
colnames(TSS_plus)[colnames(TSS_plus) == "txStart"] <- "tss"
TSS_minus <- anno[anno$strand == "-", c("name2", "chrom", "strand", "txEnd")]
colnames(TSS_minus)[colnames(TSS_minus) == "txEnd"] <- "tss"
TSSs <- rbind(TSS_plus, TSS_minus)

## Genome anno contains also non-standard chromosomes
table(TSSs$chrom)
# Select only TSSs in standard chromosomes
chroms <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")
TSSs <- TSSs[TSSs$chrom %in% chroms, ]
table(TSSs$chrom)

colnames(TSSs)[1] <- "gene_id"

#write_tsv(TSSs, "./TSSs_from_USCS_hg19_EMSEMBL.tsv")

tss <- read_tsv("/Users/ieo6983/Desktop/fragile_enhancer_clinical/data/functional_genomics/others/TSSs_elisa/TSSs_from_USCS_hg19_EMSEMBL.tsv")
length(unique((tss$gene_id)))
       