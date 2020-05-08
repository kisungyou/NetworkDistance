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
#' \donttest{
#' ## load data
#' data(graph20)
#'
#' ## compute distance matrix
#' output = nd.extremal(graph20, out.dist=FALSE, k=2)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(output$D[,20:1], main="two group case", col=gray(0:32/32), axes=FALSE)
#' par(opar)
#' }
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
