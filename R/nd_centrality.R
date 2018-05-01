#' Centrality Distance
#'
#' Centrality is a core concept in studying the topological structure of
#' complex networks, which can be either defined for each node or edge.
#' \code{nd.centrality} offers 3 distance measures on node-defined centralities.
#' See this \href{https://en.wikipedia.org/wiki/Centrality}{Wikipedia page} for more
#' on network/graph centrality.
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param mode type of node centrality definitions to be used.
#' @param directed a logical; \code{FALSE} as symmetric, undirected graph.
#'
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' \item{features}{an \eqn{(N\times M)} matrix where rows are node centralities for each graph.}
#' }
#'
#'
#' @examples
#' ## generate two types of adjacency matrices of size (3-by-3)
#' rbin1 = rbinom(9,1,0.8); mat1 = matrix(rbin1,nrow=3)
#' rbin2 = rbinom(9,1,0.2); mat2 = matrix(rbin2,nrow=3)
#'
#' mattype1 = ceiling((mat1+t(mat1))/2); diag(mattype1)=0;
#' mattype2 = ceiling((mat2+t(mat2))/2); diag(mattype2)=0;
#'
#' A = list()
#' for (i in 1:3){A[[i]]=mattype1} # first 3 are type-1
#' for (i in 4:6){A[[i]]=mattype2} # next  3 are type-2
#'
#' ## use 3 types of centrality measures
#' out1 <- nd.centrality(A,out.dist=FALSE,mode="Degree")
#' out2 <- nd.centrality(A,out.dist=FALSE,mode="Close")
#' out3 <- nd.centrality(A,out.dist=FALSE,mode="Between")
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' image(out1$D, main="Degree")
#' image(out2$D, main="Closeness")
#' image(out3$D, main="Betweenness")
#'
#' @references
#' \insertRef{roy_modeling_2014}{NetworkDistance}
#'
#' @rdname nd_centrality
#' @export
nd.centrality <- function(A, out.dist=TRUE,
                          mode=c("Degree","Close","Between"),
                          directed=FALSE){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.csd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N     = length(listA)
  M     = nrow(listA[[1]])
  # 3. out.dist & directed
  if ((!is.logical(out.dist))||(!is.logical(directed))){
    stop("* nd.centrality : 'out.dist' and 'directed' should be logical variables.")
  }
  # 4. mode
  if (missing(mode)){
    mode = "Degree"
  } else {
    mode = match.arg(mode)
  }

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. prepare for the results
  mat_features = array(0,c(N,M))
  mat_dist = array(0,c(N,N))
  #   2. transform into igraph objects & compute characteristics
  for (i in 1:N){
    #   2-1. transform
    if (directed==FALSE){
      tgt = graph_from_adjacency_matrix(listA[[i]], mode="undirected")
    } else {
      tgt = graph_from_adjacency_matrix(listA[[i]], mode="directed")
    }
    #   2-2. compute features & record
    if (mode=="Degree"){
      mat_features[i,] = as.vector(igraph::degree(tgt))
    } else if (mode=="Close"){
      mat_features[i,] = as.vector(igraph::closeness(tgt))
    } else if (mode=="Between"){
      mat_features[i,] = as.vector(igraph::betweenness(tgt))
    }
  }
  #   3. compute pairwise distances
  for (i in 1:(N-1)){
    vec1 = mat_features[i,]
    for (j in (i+1):N){
      vec2     = mat_features[j,]
      solution = sum(abs(vec1-vec2))

      mat_dist[i,j] = solution
      mat_dist[j,i] = solution
    }
  }

  #-------------------------------------------------------
  ## RETURN RESULTS
  if (out.dist){
    mat_dist = as.dist(mat_dist)
  }

  result = list()
  result$D= mat_dist
  result$features = mat_features
  return(result)
}
