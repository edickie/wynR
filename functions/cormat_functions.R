
require(igraph)

#' Calculates a correlation matrix then vectorises it
#'
#' Uses igraph to do so..
#'
#' @param data dataframe that will be correlated (column-wise)
#'
#' @return numeric vector of the flattened upper-triangle from the correlation matrix 
#'
#' @examples
#' 
#'
#' @export
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

#' Reads a mean-timeseries csv file, than calculates the roi cross-correlation and vectorises it
#'
#' Uses igraph to do so..
#'
#' @param meants_csv filepath to a csv file holding the roi timeseries data
#' @param roi.idx (optional) row.indices for the ROIs to use
#'
#' @return numeric vector of the flattened upper-triangle from the cross-correlation matrix 
#'
#' @examples
#' 
#'
#' @export
meants_to_corrZvector <- function(meants_csv, roi.idx = NULL) {
  meants <- read.csv(meants_csv, header=FALSE)
  if (!is.null(roi.idx)) {meants <- meants[roi.idx, ]}
  vec <- corrZvector(t(meants))
  return(vec)
}


#' Reads many mean-timeseries csv files, and returns a dataframe of vectorised cross-correlation values
#' 
#' Uses igraph to do so and foreach and doParallel to do it quickly
#' doParallel will use the number of cores available -2
#'
#' @param filepaths vector of filepaths to a csv file holding the roi timeseries data
#' @param myedgenames (optional) vector of edge names to use as column header for output
#' @param roi.idx (optional) row.indices for the ROIs to use (for selecting subnetworks)
#'
#' @return matrix where each row is the flattened upper-triangles from the cross-correlation matrixes 
#'
#' @examples
#' 
#'
#' @export
fc_from_meants_csvs <- function(filepaths, myedgenames=NULL, roi.idx=NULL) {
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
    meants_to_corrZvector(x, roi.idx)
  }
  if (!is.null(myedgenames)) { colnames(theZsmat) <- myedgenames}
  rownames(theZsmat) <- filepaths
  stopCluster(cl)
  return(theZsmat)
}

#' Creates a character vector of the edge names from a list of nodes 
#' 
#' Uses igraph to do so. Returns upper-triangle edge-names assuming a fully connected (or weighted) graph. 
#'
#' @param node.names character vector of node names 
#'
#' @return a character vector of edge names
#'
#' @examples
#' 
#'
#' @export
get_edgenames_uppertri <- function(node.names) {
  num.rows <- length(node.names)
  fakecormat <- matrix(rep(1, num.rows^2), 
                       nrow = num.rows, ncol = num.rows,
                       dimnames = list(node.names, node.names))
  g<-graph_from_adjacency_matrix(fakecormat,mode="upper", 
                                 weighted=T, diag=F, 
                                 add.rownames = "code")
  g.df <- as.data.frame(get.edgelist(g), names=T)
  ## get two variables of interest.. edgenames and the number of edges
  myedgenames <- paste(g.df[ ,1],g.df[ ,2],sep=".") ## the V1.V2 names
  return(myedgenames)
}

#' Determines the number of edges in the upper-triangle given the number of nodes 
#' 
#' Uses igraph to do so. Returns upper-triangle edge-names assuming a fully connected (or weighted) graph. 
#'
#' @param num.nodes integer value for the number of nodes 
#'
#' @return the number of edges in the upper-triangle
#'
#' @examples
#' 
#'
#' @export
get_length_uppertri <- function(num.nodes) {
  fakecormat <- matrix(rep(1, num.nodes^2), 
                       nrow = num.nodes, ncol = num.nodes)
  g<-graph_from_adjacency_matrix(fakecormat,mode="upper", diag=F)
  l <- length(E(g))
  return(l)
}


#' Converts a data.frame of edgewise results to matrix format 
#' 
#' Uses igraph to do so. Expects a dataframe with 3 columns: Node 1, Node2, and a vector of weights.
#'
#' @param data dataframe (3-columns) of edge names and one result column 
#'
#' @return symetrical matrix of the results
#'
#' @examples
#' 
#'
#' @export
dataframe_to_adjacency <- function(data) {
  g <- graph_from_data_frame(data, directed = F)
  gmat <- as.matrix(get.adjacency(g, type = "both"))
  return(gmat)
}