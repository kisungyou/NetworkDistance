#' Extremal distance with top-\eqn{k} eigenvalues
#'
#' Extremal distance (\code{nd.extremal}) is a type of spectral distance measures on two graphs' graph Laplacian,
#' \deqn{L := D-A}
#' where \eqn{A} is an adjacency matrix and \eqn{D_{ii}=\sum_j A_{ij}}. It takes top-\eqn{k} eigenvalues from
#' graph Laplacian matrices and take normalized sum of squared differences as metric. Note that it is
#' \emph{1. non-negative}, \emph{2. separated}, \emph{3. symmetric}, and satisfies \emph{4. triangle inequality} in that
#' it is indeed a metric.
#'
#' @param A a list of length \code{N} containing adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param k the number of largest eigenvalues to be used.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' \item{spectra}{an \eqn{(N\times k)} matrix where each row is top-\eqn{k} Laplacian eigenvalues.}
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
#' output = nd.extremal(A, out.dist=FALSE, k=2)
#' image(output$D, main="two group case")
#'
#' @references
#' \insertRef{jakobson_extremal_2002}{NetworkDistance}
#'
#' @rdname nd_extremal
#' @export
nd.extremal <- function(A, out.dist=TRUE, k=ceiling(nrow(A)/5)){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.extremal : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  # 3. k : the number of top eigenvalues to be used.
  k = as.integer(k)
  if ((length(as.vector(k))!=1)||(k<1)||(k>=nrow(listA[[1]]))){
    stop("* nd.extremal : parameter 'k' should be [1,#(nodes)).")
  }

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. setup
  N = length(listA)
  mat_eigs = array(0,c(N,k))
  mat_dist = array(0,c(N,N))

  #   2. eigenvalue computations at once
  for (i in 1:N){
    L            = as.matrix(laplacian_unnormalized(listA[[i]]))
    mat_eigs[i,] = as.vector(RSpectra::eigs(L,k)$values)
  }

  #   3. pairwise computation
  for (i in 1:(N-1)){
    eig1 = mat_eigs[i,]
    for (j in (i+1):N){
      eig2 = mat_eigs[j,]
      #   3-1. top and bottom
      numerator   = sum((eig1-eig2)^2)
      denominator = min(sum(eig1^2),sum(eig2^2))
      #   3-2. exception handling
      if (denominator==0){
        solution = NA
      } else {
        solution = sqrt(numerator/denominator)
      }
      #   3-3. plugin
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
  result$spectra = mat_eigs
  return(result)
}
