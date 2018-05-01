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
#' ## Compute Distance Matrix and Visualize
#' output = nd.hamming(A)
#' image(as.matrix(output$D), main="two group case")
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
