
library(tidyverse)
library(TFBSTools)
library(BSgenome)
#BiocManager::install("JASPAR2020")
library(JASPAR2020)
library(GenomicRanges)
library(VariantAnnotation)
#BiocManager::install(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)

SEED <- 4321
set.seed(SEED)

WIN <- 50
MARKERS <- c("CtIP", "GRHL")
motif_thresh <- 30 # score is typically expressed as a percentage of the maximum possible score for a match to the PWM. 

path_enh_SSMs <- list("CtIP" = fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_CtIP.all_overlaps.", WIN, "bp_WIN.tsv")), 
                      "GRHL" = fs::path(paste0("/Users/ieo6983/Desktop/fragile_enhancer_clinical/results/ICGC/enhancers_SSMs_overlaps/data/Table_enh_SSMs_GRHL.all_overlaps.", WIN, "bp_WIN.tsv")))


##


# GRHL2 motif ID - BaseID: MA1105; MatrixID: MA1105.2 (JASPAR2020 CORE, latest version)
motif_id <- "MA1105.2"
pfm <- getMatrixSet(JASPAR2020, opts = list(ID = motif_id))

# Load variants - enhancers overlap 
enh_SSMs <- list()
enh_gr <- list()
for(marker in MARKERS){
  enh_SSMs[[marker]] <- read_tsv(path_enh_SSMs[[marker]])
  
  # Convert 50bp extended enhancers into GRanges 
  enh_gr1 <- makeGRangesFromDataFrame(enh_SSMs[[marker]][,c(1:4,7:8)], keep.extra.columns = T, 
                                      seqnames.field = "seqnames", start.field = "start", end.field = "end", ignore.strand = T)
  enh_gr[[marker]] <- enh_gr1
}


##


## Motif scanning

# GRHL2 motif has a length of 12 bp, but matchPWM looks for complete matches (of the whole sequence)
pwm <- toPWM(pfm[[1]])

matches_list_markers <- list()
motif_pos_info <- data.frame(matrix(nrow = 1, ncol=5)) 
colnames(motif_pos_info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
motif_pos_info_markers <- list("CtIP" = motif_pos_info, "GRHL" = motif_pos_info)

for(marker in MARKERS){
  
  # Look for GRHL2 motifs matches within enhancer sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, enh_gr[[marker]])
  min_score_thresh <- paste0(as.character(motif_thresh), "%")
  matches_list <- lapply(seqs, function(seq) {
    matchPWM(pwm@profileMatrix, seq, min.score = min_score_thresh, with.score = TRUE)
  })
  matches_list_markers[[marker]] <- matches_list
  
  # Get motif start-end coordinates in enhancer regions 
  for(i in 1:length(matches_list_markers[[marker]])){
    match <- matches_list_markers[[marker]][[i]]
    if(length(match) == 0){
      # TODO: add enhancer info 
      enh_name <- enh_gr[[marker]][i]$name
      info <- data.frame(matrix(c(enh_name, rep(NA,4)), nrow = 1))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    } else {
      # CAREFUL: there can be more than 1 match for each enhancer 
      motif_start <- start(match) + start(enh_gr[[marker]][i]) - 1
      motif_end <- end(match) + end(enh_gr[[marker]][i]) - 1
      motif_score <- round(mcols(match)$score,2) 
      enh_name <- enh_gr[[marker]][i]$name
      l <- length(motif_start)
      info <- data.frame(matrix(c(rep(enh_name,l), rep(motif_id,l), motif_start, motif_end, motif_score), nrow=l))
      colnames(info) <- c("enh_name", "motif_id", "motif_start", "motif_end", "motif_score")
      motif_pos_info_markers[[marker]] <- rbind(motif_pos_info_markers[[marker]], info)
    }
  }
  
  # Print some info
  n_enh_motif <- length(unique(
    motif_pos_info_markers[[marker]] %>% 
      dplyr::filter(., ( !duplicated(.) & !is.na(motif_id) ) ) %>%
      .$enh_name
    ))
  n_enh_tot <- length(unique(enh_SSMs[[marker]]$name))
  print(paste0("Fraction of ", marker, " enhancers with a GRHL2 motif within: ", n_enh_motif, " / ", n_enh_tot, 
               " - ~", round(n_enh_motif / n_enh_tot * 100), "%"))
}


##


# TODO: 
#' This script could be used to classify mutations as damaging AFTER overlap among enhancers and SSMs
#' @Repeat motif matching with a Table_enh_SSMs overlap, i.e. with WIN = 50 bp
#' @Define 3 "damaging" levels for mutations: 
#' - 1 = within motif;
#' - 2 = outside motif, within +/- 10 bp window;
#' - 3 = outside window around motif but within +/- 50 bp from peak summit.


