
### REMOVE_SELF_OVERLAPS
# Given a GRages object containing a set of genomic regions, remove overlapping regions within the GRanges 

remove_self_overlaps <- function(input_granges){
  
  self_ovrlps <- data.frame(findOverlaps(input_granges)) # overlap GRanges_obj with itself
  self_ovrlps$hit <- self_ovrlps$queryHits == self_ovrlps$subjectHits 
  
  filter_overlaps <- self_ovrlps %>% dplyr::group_by(queryHits) %>%
    dplyr::arrange(hit, .by_group = T) %>% # arrange by T/F
    filter(row_number() == 1) %>% # and filter out regions overlapping with themselves only
    ungroup() %>% dplyr::select(hit)
  output_granges <- input_granges[filter_overlaps$hit] 
  
  return(output_granges)
}





### GENERATE RANDOM SEQS ###
# Generate a set of random genomic sequences 

generate_random_seqs <- function(seed, n_seqs, win){
  set.seed(seed)
  
  # Define chromosomes info
  chrs_info <- data.frame(
    "chrom" = c(paste("chr", seq(1:22), sep = ""), "chrX", "chrY"), 
    seq_lenghts = c(249250621,243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
                    141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392 , 90354753, 
                    81195210, 78077248, 59128983 , 63025520, 48129895, 51304566, 155270560 , 59373566)
  )
  
  # Sample chromosomes and start positions 
  ran_seqs <- data.frame(chrom = sample(chrs_info$chrom, n_seqs, replace=T))
  start_pos <- c()
  for(i in 1:length(ran_seqs$chrom)){
    chrom <- ran_seqs$chrom[i]
    chrom_length <- chrs_info[chrs_info$chrom == chrom, ]$seq_lenghts
    start_pos[i] <- round(runif(1, 0, chrom_length),0)
  }
  ran_seqs$start <- start_pos
  
  # Add end position and sequence name
  ran_seqs$end <- ran_seqs$start + (win*2)
  ran_seqs$name <- paste(ran_seqs$chrom, ran_seqs$start, ran_seqs$end, sep = "_")
  
  return(ran_seqs)  
}




### STAT ENHANCERS ###
# Computes number of donors hit by each enhancer 
# Inputs must be 2 lists: 
# - a list of data.frames with overlaps 
# - a list of the dataframes containing the data.frames containing the enhancers used for the overlaps
# classes = vector containing classes of enhancers/sequences overlapped

compute_stat_enhancers <- function(input_overlaps, source_enhancers_df, classes){
  df_all_stats <- data.frame(matrix(ncol = length(input_overlaps), 
                                    nrow = max(sapply(list(marker_enh_high, marker_enh_low), function(df) dim(df)[1]))))
  colnames(df_all_stats) <- classes
  
  for(i in 1:length(input_overlaps)){
    df_overlaps <- input_overlaps[[i]]
    df_source <- source_enhancers_df[[i]]
    
    # Compute number of donors in which each enhancer is mutated
    stat_enh <- df_overlaps %>% dplyr::group_by(name) %>%
      dplyr::summarise(., n_donors_hit = length(unique(icgc_donor_id)))
    
    # Add to df enhancers with no mutation in any patient
    stat_enh <- rbind(stat_enh, 
                      data.frame("name" = df_source$name[!df_source$name %in% stat_enh$name], 
                                 "n_donors_hit" = rep(0, length(df_source$name[!df_source$name %in% stat_enh$name]))))
    
    # add to df_all_stats
    df_all_stats[, i] <- c(stat_enh$n_donors_hit, 
                           rep(NA, length(df_all_stats[, i]) - length(stat_enh$n_donors_hit)))
  }
  
  return(df_all_stats)
}





### STAT DONORS
# Compute number of enhancers hit/mutated in each donor
# input_overlaps = df of overlaps
# source_donors_df = df with all donors used up to now

compute_stat_donors <- function(input_overlaps, source_donors_df, classes){
  df_all_stats <- data.frame(matrix(ncol = length(input_overlaps), 
                                    nrow = length(unique(source_donors_df$icgc_donor_id))))
  colnames(df_all_stats) <- classes
  
  for(i in 1:length(input_overlaps)){
    df_overlaps <- input_overlaps[[i]]
    
    # Compute number of enhancers mutated in each donor
    stat_donors <- df_overlaps %>% group_by(icgc_donor_id) %>%
      dplyr::summarise("n_enhancers_hit" = length(unique(name)))
    
    # Add to df donors in which no enhancer is mutated
    stat_donors <- rbind(stat_donors, 
                         data.frame("icgc_donor_id" = unique(source_donors_df$icgc_donor_id[!source_donors_df$icgc_donor_id %in% stat_donors$icgc_donor_id]), 
                                    "n_enhancers_hit" = rep(0, length(unique(source_donors_df$icgc_donor_id[!source_donors_df$icgc_donor_id %in% stat_donors$icgc_donor_id])))))
    # add to df_all_stats
    df_all_stats[, i] <- stat_donors$n_enhancers_hit
  }
  
  return(df_all_stats)
}

