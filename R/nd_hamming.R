#' Hamming Distance
#'
#' Hamming Distance is the count of discrepancy between two binary networks for each edge.
#' Therefore, if used with non-binary networks, it might return a warning message and distorted results.
#' It was originally designed to compare two strings of equal length, see \href{https://en.wikipedia.org/wiki/Hamming_distance}{Wikipedia page} for more detailed introduction.
#'
#' @param A a list of length \eqn{N} containing adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' }
#'
#' @examples
#' ## load example data and extract only a few
#' data(graph20)
#' gr.small = graph20[c(1:5,11:15)]
#'
#' ## compute distance matrix
#' output = nd.hamming(gr.small, out.dist=FALSE)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(output$D[,10:1], main="two group case", axes=FALSE, col=gray(0:32/32))
#' par(opar)
#'
#' @references
#' \insertRef{hamming_error_1950}{NetworkDistance}
#'
#' @rdname nd_hamming
#' @export
nd.hamming <- function(A, out.dist=TRUE){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.csd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. setup
  N = length(listA)
  M = nrow(listA[[1]])
  mat_dist = array(0,c(N,N))

  #   2. warning
  for (i in 1:N){
    if (length(unique(as.vector(listA[[i]])))>2){
      message("* nd.hamming : Hamming distance is originally used for binary networks. The result may not be valid.")
      break
    }
  }

  #   3. pairwise distance
  for (i in 1:(N-1)){
    mat1 = listA[[i]]
    diag(mat1)=0
    for (j in (i+1):N){
      mat2 = listA[[j]]
      diag(mat2)=0

      solution = sum(abs(mat1-mat2))/(M*(M-1))
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
  return(result)
}
