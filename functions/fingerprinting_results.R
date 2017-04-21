#' Calculates "fingerprinting" style stats from a subject to subject correlation matrix 
#' 
#' Uses igraph to do so.  
#'
#' @param sub2sub_cormat correlation matrix of M X N subjects with subject ID as the row and column names 
#'
#' @return dataframe of within and between subject values. 
#' 
#' Each row representing stats calculated from one subject's row of the input matrix.
#' Results include.
#'    subid1 the subject ID, taken from the row.name of the input matrix 
#'    within_popZ   The Z score for the within subject correlation within the row.
#'    within_rank  The percentile rank of the withing subject correlation within the row.
#'    within_Zrho      The within subject correlation (Z-transformed)
#'    cross_meanZ       The mean of the cross subject (Z-transformed) correlation (within that row)
#'    cross_sdZ         The stanard deviation of the cross subject (Z-transformed) correlation (within that row)
#' 
#' @details 
#' Currently only works when there is one subject ID in the row.names that matches one subject ID in teh column names.
#' Requires dplyr and tidyr
#'
#' @export
fingerprint_results <- function(sub2sub_cormat) {
  require(dplyr)
  require(tidyr)
  tmpdf <- as.data.frame(sub2sub_cormat) # convert to dataframe
  tmpdf$subid1 <- row.names(tmpdf)       # make subid1 column from row.names
  # melt the datafrom, Z-transform the correlations, 
  # and label within and cross connections based on if original row and column names match
  result <- tmpdf %>%
    gather(subid2, rho, -subid1) %>%    
    mutate(Z = atanh(rho),
           within_cross = if_else(subid1==subid2, "within", "cross")) %>%
    ## calculate population Z score and rank within each row
    group_by(subid1) %>%
    mutate(popZ = as.numeric(scale(Z)),
           rankZ = as.numeric(rank(Z), ties.method = "max")/n()) %>%
    ## calculate summary statistics for both within and cross participant values
    ungroup() %>%
    group_by(subid1,within_cross) %>%
    summarise(n = n(),
              meanZ = mean(Z),
              sdZ = sd(Z),
              med_popZ = median(popZ),
              med_rankZ = median(rankZ)) %>%
    ungroup() %>%
    ## rearrange the output to a one row per participant format
    gather(measure, value, -subid1, -within_cross) %>%
    unite(mycolnames, within_cross, measure, sep = "_") %>%
    spread(mycolnames, value) %>%
    ## remove participants where the their was not matching column name from the result
    filter(within_n == 1, cross_n == (ncol(sub2sub_cormat)-1)) %>%
    select(subid1, 
           within_med_popZ, within_med_rankZ, within_meanZ, 
           cross_meanZ, cross_sdZ)
  ## rename some of the columns
  names(result) <- c("subid", "within_popZ", "within_rank", "within_Zrho", 
                     "cross_meanZ", "cross_sdZ")
  return(result)
}