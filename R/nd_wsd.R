#' Distance with Weighted Spectral Distribution
#'
#' Normalized Laplacian matrix contains topological information of
#' a corresponding network via its spectrum. \code{nd.wsd} adopts weighted
#' spectral distribution of eigenvalues and brings about a metric via
#' binning strategy.
#'
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param K the number of bins for the spectrum interval \eqn{[0,2].}
#' @param wN a decaying exponent; default is \eqn{4} set by authors.
#'
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' \item{spectra}{an \eqn{(N\times M)} matrix of rows being eigenvalues for each graph.}
#' }
#'
#'
#' @examples
#' ## generate two types of adjacency matrices of size (3-by-3)
#' rbin1 = rbinom(9,1,0.8); mat1 = matrix(rbin1,nrow=3)
#' rbin2 = rbinom(9,1,0.2); mat2 = matrix(rbin2,nrow=3)
#'
#' mattype1 = ceiling((mat1+t(mat1))/2)
#' mattype2 = ceiling((mat2+t(mat2))/2)
#'
#' A = list()
#' for (i in 1:3){A[[i]]=mattype1} # first 3 are type-1
#' for (i in 4:6){A[[i]]=mattype2} # next  3 are type-2
#'
#' ## Compute Distance Matrix and Visualize
#' output = nd.wsd(A, out.dist=FALSE, K=10)
#' image(output$D, main="two group case")
#'
#' @references
#' \insertRef{fay_weighted_2010}{NetworkDistance}
#'
#' @rdname nd_wsd
#' @export
nd.wsd <- function(A, out.dist=TRUE, K=50, wN=4){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.wsd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N     = length(listA)    # number of networks
  M     = nrow(listA[[1]]) # number of nodes
  # 3. K : the number of binning
  K = as.integer(K)
  if ((length(as.vector(K))!=1)||(K<1)||(is.na(K))||(is.infinite(K))){
    stop("* nd.wsd : parameter 'K' for binning should be a positive integer.")
  }
  # 4. wN : weight function exponent
  wN = as.integer(wN)
  if ((length(as.vector(wN))!=1)||(wN<1)||(is.na(wN))||(is.infinite(wN))){
    stop("* nd.wsd : weight exponent 'wN' is a positive integer; 4 was chosen by authors.")
  }

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. setup
  mat_eigs = array(0,c(N,M))
  mat_dist = array(0,c(N,N))

  #   2. eigenvalue computation
  for (i in 1:N){
    mat_eigs[i,] = as.vector(eigen(as.matrix(laplacian_normalized(listA[[i]])))$values)
  }

  #   3. binning and computing the distance
  testgrid  = seq(from=0,to=2,length.out = (K+1))
  middlepts = testgrid[1:K]+(diff(testgrid)/2)
  for (i in 1:(N-1)){
    spect1 = mat_eigs[i,]
    bhist1 = wsd_binhist(spect1, testgrid)
    for (j in (i+1):N){
      spect2   = mat_eigs[j,]
      bhist2   = wsd_binhist(spect2, testgrid)

      # compute distance given two counting histograms and middle points for bins
      solution = wsd_histdist(bhist1, bhist2, middlepts, wN)

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



#  ------------------------------------------------------------------------
#  generate counting histogram given grids
#' @keywords internal
#' @noRd
wsd_binhist <- function(eig, grid){
  K = length(grid)-1
  count = rep(0,K)
  for (i in 1:K){
    count[i] = length(intersect(which(eig>=grid[i]),which(eig<grid[i+1])))
  }
  count = count/sum(count)
  return(count)
}

# given two histograms and middle points, compute distance as D5.
#' @keywords internal
#' @noRd
wsd_histdist <- function(h1, h2, mpts, wN){
  return(sum(((1-mpts)^wN)*((h1-h2)^2)))
}




