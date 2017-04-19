
require(igraph)

corrZvector <- function(data) {
  ## correlate and graph
  cormat <- cor(data)
  g<-graph_from_adjacency_matrix(cormat,mode="upper", 
                                 weighted=T, diag=F)
  # take the egde list as a vector
  thecorrs <- E(g)$weight
  
  # apply the Z transform (so we can do stats)
  theseZ <- atanh(thecorrs)
  
  # save the output to the data.frame
  return(theseZ)
}

meants_to_corrZvector <- function(meants_csv) {
  meants <- read.csv(meants_csv, header=FALSE)
  vec <- corrZvector(t(meants))
  return(vec)
}

fc_from_meants_csvs <- function(filepaths, myedgenames) {
  require(foreach)
  library(doParallel)
  no_cores <- detectCores()
  if (no_cores > 3) {no_cores <- no_cores - 2}
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl)
  theZsmat <- foreach(x=filepaths, 
                      .combine = rbind,
                      .export=c('meants_to_corrZvector', 'corrZvector'), 
                      .packages='igraph'
                      ) %dopar% {
    meants_to_corrZvector(x)
  }
  colnames(theZsmat) <- myedgenames
  rownames(theZsmat) <- filepaths
  stopCluster(cl)
  return(theZsmat)
}


get_edgenames_uppertri <- function(row.names) {
  num.rows <- length(row.names)
  fakecormat <- matrix(rep(1, num.rows^2), 
                       nrow = num.rows, ncol = num.rows,
                       dimnames = list(row.names, row.names))
  g<-graph_from_adjacency_matrix(fakecormat,mode="upper", 
                                 weighted=T, diag=F, 
                                 add.rownames = "code")
  g.df <- as.data.frame(get.edgelist(g), names=T)
  ## get two variables of interest.. edgenames and the number of edges
  myedgenames <- paste(g.df[ ,1],g.df[ ,2],sep=".") ## the V1.V2 names
  return(myedgenames)
}


fingerprint_results <- function(sub2sub_cormat) {
  require(dplyr)
  require(tidyr)
  tmpdf <- as.data.frame(sub2sub_cormat)
  tmpdf$subid1 <- row.names(tmpdf)
  result <- tmpdf %>%
    gather(subid2, rho, -subid1) %>%
    mutate(Z = atanh(rho),
           within_cross = if_else(subid1==subid2, "within", "cross")) %>%
    group_by(subid1) %>%
    mutate(popZ = as.numeric(scale(Z)),
          rankZ = as.numeric(rank(Z), ties.method = "first")/n()) %>%
    ungroup() %>%
    group_by(subid1,within_cross) %>%
    summarise(n = n(),
              meanZ = mean(Z),
              sdZ = sd(Z),
              med_popZ = median(popZ),
              med_rankZ = median(rankZ)) %>%
    ungroup() %>%
    gather(measure, value, -subid1, -within_cross) %>%
    unite(mycolnames, within_cross, measure, sep = "_") %>%
    spread(mycolnames, value) %>%
    filter(within_n == 1, cross_n == (ncol(sub2sub_cormat)-1)) %>%
    select(subid1, 
           within_med_popZ, within_med_rankZ, within_meanZ, 
           cross_meanZ, cross_sdZ)
  return(result)
}

dataframe_to_adjacency <- function(data) {
  g <- graph_from_data_frame(newtable, directed = F)
  gmat <- as.matrix(get.adjacency(g, type = "both"))
  return(gmat)
}