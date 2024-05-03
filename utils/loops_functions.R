
### --- Functions --- ###


#' Find and filter loops overlapping with enhancers - Method 1 
#'
#' This function annotates loops' bins to a set of input enhancers
#' 
#' @param enhancer_file File with enhancers coordinates: chrom, start, end, summit, name, (cluster). I.e. GRHL2-enhancers
#' @param loops_file File with loops coordinates: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', others (p-value, FDR, ..)
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param ext_bins Number of bins left/right by which loop bin is extended. Overlaps are searched within the extended bin.
#' @returns Input loops_file with 4 additional columns: name1, name2, cluster1, cluster2. 
#' @examples
#' scr_enh_overlaps <- find_loops_overlapping_enhancers(enhancers_file = enh_grhl2, loops_file = scr, bin_len_kb = 4, ext_nbins = 1)

find_loops_overlapping_enhancers <- function(enhancers_file, 
                                             loops_file, 
                                             bin_len_kb = 4,
                                             ext_nbins = 1){
  
  # Conver input files to GRanges objects
  enhancers_gr <- makeGRangesFromDataFrame(enhancers_file, keep.extra.columns = T)
  loops_bin1 <- loops_file[, endsWith(colnames(loops_file), "1")]
  loops_bin2 <- loops_file[, endsWith(colnames(loops_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  # Check if all regions have equal width, equal to kb
  if(all(width(loops_bin1_gr) == width(loops_bin1_gr)[1])){
    message(sprintf("All loop bins 1 have equal length; equal to: %d", width(loops_bin1_gr)[1]))
  } else {
    message("WARNING: not all loop bins 1 have equal length")
  }
  
  if(all(width(loops_bin2_gr) == width(loops_bin2_gr)[1])){
    message(sprintf("All loop bins 2 have equal length; equal to: %d", width(loops_bin2_gr)[1]))
  } else {
    message("WARNING: not all loop bins 1 have equal length")
  }
  
  # Extend loop bins of +/- n bins
  start(loops_bin1_gr) <- start(loops_bin1_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin1_gr) <- end(loops_bin1_gr) + (bin_len_kb*1000*ext_nbins)
  start(loops_bin2_gr) <- start(loops_bin2_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin2_gr) <- end(loops_bin2_gr) + (bin_len_kb*1000*ext_nbins)
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, enhancers_gr, select = "arbitrary")
  enh_1 <- enhancers_file[ovrlp_bin1, c("name", "cluster")]
  colnames(enh_1) <- paste0(colnames(enh_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, enh_1) # GRanges are extended, but df remains as originally
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, enhancers_gr, select = "arbitrary")
  enh_2 <- enhancers_file[ovrlp_bin2, c("name", "cluster")]
  colnames(enh_2) <- paste0(colnames(enh_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, enh_2) # GRanges are extended, but df remains as originally
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>%
    relocate(., c(name1,cluster1), .after = end2)
  
  return(loops_ovrlp)
}


#' Find and filter loops overlapping with enhancers - Method 2
#'
#' This function annotates loops' bins to a set of input enhancers
#' 
#' @param enhancer_file File with enhancers coordinates: chrom, start, end, summit, name, (cluster). I.e. GRHL2-enhancers
#' @param loops_file File with loops coordinates: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', others (p-value, FDR, ..)
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param kb_win  Window in kb to extend enhancers (from summit). Overlaps with enhancers are searched considering the extended window.
#' @returns Input loops_file with 4 additional columns: name1, name2, cluster1, cluster2. 
#' @examples
#' scr_enh_overlaps <- find_loops_overlapping_enhancers_m2(enhancers_file = enh_grhl2, loops_file = scr, bin_len_kb = 4, ext_nbins = 1)

find_loops_overlapping_enhancers_m2 <- function(enhancers_file, 
                                             loops_file, 
                                             bin_len_kb = 4,
                                             kb_win = 1){
  
  # Extend enhancer regions and conver to GRanges
  enhancers_file$start <- enhancers_file$summit - (kb_win * 1000)
  enhancers_file$end <- enhancers_file$summit + (kb_win * 1000)
  enhancers_gr <- makeGRangesFromDataFrame(enhancers_file, keep.extra.columns = T)
  
  # Extract bins and convert to GRanges
  loops_bin1 <- loops_file[, endsWith(colnames(loops_file), "1")]
  loops_bin2 <- loops_file[, endsWith(colnames(loops_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, enhancers_gr, select = "arbitrary")
  enh_1 <- enhancers_file[ovrlp_bin1, c("name", "cluster")]
  colnames(enh_1) <- paste0(colnames(enh_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, enh_1) 
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, enhancers_gr, select = "arbitrary")
  enh_2 <- enhancers_file[ovrlp_bin2, c("name", "cluster")]
  colnames(enh_2) <- paste0(colnames(enh_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, enh_2) 
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>%
    relocate(., c(name1,cluster1), .after = end2)
  
}


#' Find and filter loops overlapping with gene promoters, deprecated
#' 
#' Given loops already annotated to enhancers, this function annotates loops' bins to a set of input gene promoters 
#' Deprecated, use find_additional_DEGS_overlap instead
#' 
#' @param tss_file File with TSSs coordinates: gene_id, chrom, strand, tss
#' @param loops_enh_file File with coordinates of loops, with info on overlapping enhancers: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', name1, cluster1, ...
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param ext_nbins Number of bins left/right by which loop bin is extended. Overlaps are searched within the extended bin.
#' @returns Input loops_enh_file with 8 additional columns: gene_name1, DE1, log2FoldChange1, padj1, (gene_name2, DE2, log2FoldChange2, padj2)

find_additional_promoters_overlap <- function(tss_file, 
                                             loops_enh_file, 
                                             bin_len_kb = 4,
                                             ext_nbins = 1){
  
  # Create GRanges obj for TSSs
  tss_file$start <- tss_file$tss
  tss_file$end <- tss_file$tss
  tss_gr <- makeGRangesFromDataFrame(tss_file, keep.extra.columns = T)
  
  # Create GRanges obj for loops_enh
  loops_bin1 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "1")]
  loops_bin2 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  
  # Extend loop bins of +/- n bins
  start(loops_bin1_gr) <- start(loops_bin1_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin1_gr) <- end(loops_bin1_gr) + (bin_len_kb*1000*ext_nbins)
  start(loops_bin2_gr) <- start(loops_bin2_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin2_gr) <- end(loops_bin2_gr) + (bin_len_kb*1000*ext_nbins)
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, tss_gr, select = "arbitrary")
  tss_1 <- tss_file[ovrlp_bin1, c("gene_name")]
  colnames(tss_1) <- paste0(colnames(tss_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, tss_1) # GRanges are extended, but df remains as originally
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, tss_gr, select = "arbitrary")
  tss_2 <- tss_file[ovrlp_bin2, c("gene_name")]
  colnames(tss_2) <- paste0(colnames(tss_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, tss_2) # GRanges are extended, but df remains as originally
  
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>%
    relocate(., c(name1,cluster1,gene_name1), .after = end2)
  
  return(loops_ovrlp)
}


##


#' Find and filter loops overlapping with DEGs promoters - Method 1 
#' 
#' Given loops already annotated to enhancers, this function annotates loops' bins to a set of input DEGs promoters 
#' 
#' @param DEGs_tss_file File with TSSs coordinates: gene_id, chrom, strand, tss
#' @param loops_enh_file File with coordinates of loops, with info on overlapping enhancers: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', name1, cluster1, ...
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param ext_nbins Number of bins left/right by which loop bin is extended. Overlaps are searched within the extended bin.
#' @returns Input loops_enh_file with 8 additional columns: gene_name1, DE1, log2FoldChange1, padj1, (gene_name2, DE2, log2FoldChange2, padj2)

find_additional_DEGS_overlap <- function(DEGs_tss_file, 
                                         loops_enh_file, 
                                         bin_len_kb = 4,
                                         ext_nbins = 1){
  
  tss_file <- DEGs_tss_file
  # Create GRanges obj for TSSs
  tss_file$start <- tss_file$tss
  tss_file$end <- tss_file$tss
  tss_file <- tss_file %>% dplyr::select(gene_name, chrom, tss, start, end, DE, log2FoldChange, padj)
  tss_gr <- makeGRangesFromDataFrame(tss_file, keep.extra.columns = T, seqnames.field = "chrom")
  
  # Create GRanges obj for loops_enh
  loops_bin1 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "1")]
  loops_bin2 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  
  # Extend loop bins of +/- n bins
  start(loops_bin1_gr) <- start(loops_bin1_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin1_gr) <- end(loops_bin1_gr) + (bin_len_kb*1000*ext_nbins)
  start(loops_bin2_gr) <- start(loops_bin2_gr) - (bin_len_kb*1000*ext_nbins)
  end(loops_bin2_gr) <- end(loops_bin2_gr) + (bin_len_kb*1000*ext_nbins)
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, tss_gr, select = "arbitrary")
  tss_1 <- tss_file[ovrlp_bin1, c("gene_name", "DE", "log2FoldChange", "padj")]
  colnames(tss_1) <- paste0(colnames(tss_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, tss_1) # GRanges are extended, but df remains as originally
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, tss_gr, select = "arbitrary")
  tss_2 <- tss_file[ovrlp_bin2, c("gene_name", "DE", "log2FoldChange", "padj")]
  colnames(tss_2) <- paste0(colnames(tss_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, tss_2) # GRanges are extended, but df remains as originally
  
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>% 
    relocate(., c(name1,cluster1,gene_name1, DE1, log2FoldChange1, padj1), .after = end2)
  # Remove eventual overlaps
  message(sprintf("Removing %d duplicated loops", sum(duplicated(loops_ovrlp))))
  loops_ovrlp <- loops_ovrlp[!duplicated(loops_ovrlp), ]
  
  return(loops_ovrlp)
}


##


#' Find and filter loops overlapping with DEGs promoters - Method 2
#' 
#' Given loops already annotated to enhancers, this function annotates loops' bins to a set of input DEGs promoters 
#' 
#' @param DEGs_tss_file File with TSSs coordinates: gene_id, chrom, strand, tss
#' @param loops_enh_file File with coordinates of loops, with info on overlapping enhancers: 'seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2', name1, cluster1, ...
#' @param bin_len_kb Length of loop bins resulting from loop calling in kb. I.e. 4, for 4 kb.
#' @param kb_win Window in kb to extend DEGs promoters (from TSS). Overlaps with DEGs are searched considering the extended window
#' @returns Input loops_enh_file with 8 additional columns: gene_name1, DE1, log2FoldChange1, padj1, (gene_name2, DE2, log2FoldChange2, padj2)

find_additional_DEGS_overlap_m2 <- function(DEGs_tss_file, 
                                            loops_enh_file, 
                                            bin_len_kb = 4,
                                            kb_win = 1){
  
  # Extend TSSs regions and conver to GRanges
  tss_file <- DEGs_tss_file
  tss_file$start <- tss_file$tss - (kb_win*1000)
  tss_file$end <- tss_file$tss + (kb_win*1000)
  tss_file <- tss_file %>% dplyr::select(gene_name, chrom, tss, start, end, DE, log2FoldChange, padj)
  tss_gr <- makeGRangesFromDataFrame(tss_file, keep.extra.columns = T, seqnames.field = "chrom")
  
  # Create GRanges obj for loops_enh
  loops_bin1 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "1")]
  loops_bin2 <- loops_enh_file[, endsWith(colnames(loops_enh_file), "2")]
  cols <- colnames(loops_bin1)
  loops_bin1_gr <- makeGRangesFromDataFrame(loops_bin1, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  cols <- colnames(loops_bin2)
  loops_bin2_gr <- makeGRangesFromDataFrame(loops_bin2, seqnames.field = cols[1], 
                                            start.field = cols[2], end.field = cols[3])
  
  
  # Find overlaps 
  #bin1
  ovrlp_bin1 <- findOverlaps(loops_bin1_gr, tss_gr, select = "arbitrary")
  tss_1 <- tss_file[ovrlp_bin1, c("gene_name", "DE", "log2FoldChange", "padj")]
  colnames(tss_1) <- paste0(colnames(tss_1), 1)
  loops_bin1_ovrlp <- cbind(loops_bin1, tss_1) 
  #bin2
  ovrlp_bin2 <- findOverlaps(loops_bin2_gr, tss_gr, select = "arbitrary")
  tss_2 <- tss_file[ovrlp_bin2, c("gene_name", "DE", "log2FoldChange", "padj")]
  colnames(tss_2) <- paste0(colnames(tss_2), 2)
  loops_bin2_ovrlp <- cbind(loops_bin2, tss_2)
  
  #joined
  loops_ovrlp <- cbind(loops_bin1_ovrlp, loops_bin2_ovrlp) %>% 
    relocate(., c(name1,cluster1,gene_name1, DE1, log2FoldChange1, padj1), .after = end2)
  # Remove eventual overlaps
  message(sprintf("Removing %d duplicated loops", sum(duplicated(loops_ovrlp))))
  loops_ovrlp <- loops_ovrlp[!duplicated(loops_ovrlp), ]
  
  return(loops_ovrlp)
}


##


#' Filter out ambiguous loops, deprecated
#' 
#' This function allows to modify the annotation of loops annotated to multiple regions (ambiguous) to a non-ambiguous annotation
#' Deprecated, use unbundle_ambiguous_bins instead
#' 
#' @param df_loops_enh_prom File with loops annotated to both enhancers and DEGs promoters - with ambiguously annotated loops
#' @returns Input enh_degs_loops_file with unique ENH-DEG annotation. The annotation of ambiguous is changed. 
 
filter_ambiguous_loops <- function(df_loops_enh_prom){
  
  # Annotation is ambiguous if either bin1 or bin2 overlap with both enh and prom
  cond_root_amb <- !is.na(df_loops_enh_prom$name1) & !is.na(df_loops_enh_prom$gene_name1) | !is.na(df_loops_enh_prom$name2) & !is.na(df_loops_enh_prom$gene_name2)
  
  # Split loops_enh_prom into ambiguous and not_ambiguous loops
  df_not_amb <- df_loops_enh_prom[!cond_root_amb, ]
  df_amb <- df_loops_enh_prom[cond_root_amb, ]
  message(sprintf("Check if size of df_amb and df_not_amb correspond to size of whole: %s", 
          dim(df_not_amb)[1] + dim(df_amb)[1] == dim(df_loops_enh_prom)[1]))
  
  # Store whether bin1 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin1 <- rowSums(!is.na(df_amb %>% dplyr::select(name1, gene_name1)))
  # Store whether bin2 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin2 <- rowSums(!is.na(df_amb %>% dplyr::select(name2, gene_name2)))
  
  ## Case 2-0: If bin1 overlaps with enhancer and promoter, and bin2 with "any", keep enhancer
  df_amb[df_amb$n_bin1 == 2 & df_amb$n_bin2 == 0, ]$gene_name1 <- NA
  df_amb[df_amb$n_bin1 == 0 & df_amb$n_bin2 == 2, ]$gene_name2 <- NA
  
  ## Case 2-1: If bin1 overlaps with enh and prom, and bin2 with enh|prom. In bin1, keep opposite region wrt bin2.
  # Option 1:
  cond_case2 <- df_amb$n_bin1 == 2 & df_amb$n_bin2 == 1
  # bin2 overlaps with enhancer - bin1 keep promoter, remove enhancers
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$name2),]$name1 <- NA
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$name2),]$cluster1 <- NA
  # bin2 overlaps with promoter - bin1 keep enhancer, remove promoter
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$gene_name2),]$gene_name1 <- NA
  
  # Option 2:
  cond_case2 <- df_amb$n_bin1 == 1 & df_amb$n_bin2 == 2
  # bin1 overlaps with enhancer - bin2 keep promoter, remove enhancers
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$name1),]$name2 <- NA
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$name1),]$cluster2 <- NA
  # bin1 overlaps with promoter - bin2 keep enhancer, remove promoter
  df_amb[cond_case2, ][!is.na(df_amb[cond_case2, ]$gene_name1),]$gene_name2 <- NA
  
  ## Case 2-2: if both bin overlap with enh&prom. 
  #If promoters and enhancers are all equal, bin1 keep enh, bin2 keep prom.
  #If promoters and enhancer are not equal, drop loops and output number of droppings.
  cond_case3 <- df_amb$n_bin1 == 2 & df_amb$n_bin2 == 2
  # All equal
  cond_case3_2 <- df_amb$name1 == df_amb$name2 & df_amb$gene_name1 == df_amb$gene_name2 
  df_amb[cond_case3 & cond_case3_2, ]$gene_name1 <- NA # Keep only enh in bin1
  df_amb[cond_case3 & cond_case3_2, ]$name2 <- NA # Keep only prom in bin2
  df_amb[cond_case3 & cond_case3_2, ]$cluster2 <- NA
  # Not all equal
  n_drop <- dim(df_amb[(cond_case3 & !cond_case3_2),])[1]
  message(sprintf("Dropping %d over %d loops due to ambiguous annotation", n_drop, dim(df_amb)[1]))
  df_amb <- df_amb[!(cond_case3 & !cond_case3_2),] 
  
  # In the end all loops should be 1-1 or 1-0
  df_amb$n_bin1 <- rowSums(!is.na(df_amb %>% dplyr::select(name1, gene_name1)))
  df_amb$n_bin2 <- rowSums(!is.na(df_amb %>% dplyr::select(name2, gene_name2)))
  message(sprintf("After filtering, all bins have unique annotation? %s", sum(df_amb$n_bin1 == 2 | df_amb$n_bin2 == 2) == 0))
  
  # Re-join df_not_amb with df_amb filtered
  df_amb$n_bin1 <- NULL
  df_amb$n_bin2 <- NULL
  df_loops_enh_prom <- rbind(df_not_amb, df_amb)
  
  return(df_loops_enh_prom)
}


##


#' Unbundle ambigously annotated bins 
#' 
#' This function allows to split loops annotated to multiple regions (ambiguous) into multiple non-ambiguous loops
#' 
#' @param enh_degs_loops_file File with loops annotated to both enhancers and DEGs promoters - with ambiguously annotated loops
#' @returns Input enh_degs_loops_file with unique ENH-DEG annotation. Ambiguous loops are unbundled, resulting in additional rows.
#' Ambiguous loop: ENH1/DEG1 – ENH2/DEG2 --> Split it into 2 separate loops: (1) ENH1 – DEG2, (2) DEG1 – ENH2. 

unbundle_ambiguous_bins <- function(enh_degs_loops_file){
  
  # Count how many overlapping regions x pair of bins
  n <- rowSums(!is.na(enh_degs_loops_file[, c("name1", "gene_name1", "name2", "gene_name2")]))
  cond_amb <- n>2
  
  df_not_amb <- enh_degs_loops_file[!cond_amb, ]
  df_amb <- enh_degs_loops_file[cond_amb, ]
  message(sprintf("Check if size of df_amb and df_not_amb correspond to size of whole: %s", 
                  dim(df_not_amb)[1] + dim(df_amb)[1] == dim(enh_degs_loops_file)[1]))
  
  # Store whether bin1 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin1 <- rowSums(!is.na(df_amb %>% dplyr::select(name1, gene_name1)))
  # Store whether bin2 overlaps with 1 or 2 regions (enhancer & promoter)
  df_amb$n_bin2 <- rowSums(!is.na(df_amb %>% dplyr::select(name2, gene_name2)))
  
  # Case 1-2: bin2 overlaps with both prom and enh - Case 2-1: bin1 overlaps with both prom and enh
  df_amb_new <- df_amb
  
  for(i in 1:dim(df_amb)[1]){
    r <- df_amb[i, ]
    
    # If case 1, duplicate row separating prom and enh in bin2
    if(r$n_bin1 == 1 & r$n_bin2 == 2){
      df_amb_new[i, ]$gene_name2 <- NA
      df_amb_new[i, c("DE2", "log2FoldChange2", "padj2")] <- NA
      
      r$name2 <- NA
      r$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r)
    }
    
    if(r$n_bin1 == 2 & r$n_bin2 == 1){
      df_amb_new[i, ]$gene_name1 <- NA
      df_amb_new[i, c("DE1", "log2FoldChange1", "padj1")] <- NA
      
      r$name1 <- NA
      r$cluster1 <- NA
      df_amb_new <- rbind(df_amb_new, r)
    }
    
    if(r$n_bin1 == 2 & r$n_bin2 == 2){
      r1 <- r; r2 <- r; r3 <- r
      
      df_amb_new[i, ]$gene_name1 <- NA; df_amb_new[i, ]$gene_name2 <- NA
      df_amb_new[i, c("DE1", "log2FoldChange1", "padj1")] <- NA; df_amb_new[i, c("DE2", "log2FoldChange2", "padj2")] <- NA
      
      r1$gene_name1 <- NA; r1$name2 <- NA
      r1[, c("DE1", "log2FoldChange1", "padj1")] <- NA
      r1$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r1)
      
      r2$name1 <- NA; r2$gene_name2 <- NA
      r2$cluster1 <- NA
      r2[, c("DE2", "log2FoldChange2", "padj2")] <- NA
      df_amb_new <- rbind(df_amb_new, r2)
      
      r3$name1 <- NA; r3$name2 <- NA
      r3$cluster1 <- NA; r3$cluster2 <- NA
      df_amb_new <- rbind(df_amb_new, r3)
      
    }
  }
  
  # Check if size of df_amb_new matches df_amb with duplicated rows
  check <- dim(df_amb)[1]*2
  check <- check + (dim(df_amb%>%filter(n_bin1==2&n_bin2==2))[1]*2)
  message(sprintf("Size of new df is compatible with input: %s", check == dim(df_amb_new)[1]))
  
  # Make sure that now no ambiguous bins are present
  df_amb_new$n_bin1 <- rowSums(!is.na(df_amb_new %>% dplyr::select(name1, gene_name1)))
  df_amb_new$n_bin2 <- rowSums(!is.na(df_amb_new %>% dplyr::select(name2, gene_name2)))
  message(sprintf("No more ambiguous bins: %s", sum(c(df_amb_new$n_bin1, df_amb_new$n_bin2) > 1)==0))
  
  # Keep only enh-degs loops - Remove from df_amb news rows representing loops enh-any or prom-any 
  df_amb_new_filt <- df_amb_new %>% filter((!is.na(name1) & !is.na(gene_name2)) | (!is.na(gene_name1) & !is.na(name2)))
  df_amb_new_filt$n_bin1 <- NULL; df_amb_new_filt$n_bin2 <- NULL
  
  # Re-join df_not_amb & df_amb
  df_new <- rbind(df_not_amb, df_amb_new_filt)
  
  return(df_new)
}


#' GSEA Analysis
#'
#' This function performs a GSEA analysis 
#' 
#' @param gene_list Ranked gene list ( numeric vector, names of vector should be gene names). I.e. ranked by log2FC or padj.
#' @param GO_file File representing one gene set from the Human Molecular Signatures Database (.gmt)
#' @param p_val P-value threshold used to filter enriched pathways
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @returns this returns a data.frame and a plot with enriched pathways 
#' @examples
#' hallmark <- "GSEA/DB/hallmark_gene_sets.h.all.v2023.2.Hs.symbols.gmt"
#' ranked_genes_df <- df_LFC_sig %>% arrange(., -log2FoldChange)
#' ranked_genes <- ranked_DEGs_df$log2FoldChange; names(ranked_genes) <- ranked_genes_df$gene_name;
#' GSEA(gene_list = ranked_genes, GO_file = hallmark, p_val = 0.05)

GSEA = function(gene_list, GO_file, p_val=0.05, minSize=15, maxSize=500, collapse = T) {
  set.seed(4321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=minSize, 
                        maxSize=maxSize) %>% 
    dplyr::filter(pval <= p_val) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  if(collapse == T){
    message("Collapsing Pathways -----")
    concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                        pathways = myGO,
                                        stats = gene_list)
    fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
    message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  }
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header)
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

