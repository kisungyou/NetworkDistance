#' Discrete Spectral Distance
#'
#' Discrete Spectral Distance (DSD) is defined as the Euclidean distance between
#' the spectra of various matrices, such as adjacency matrix \eqn{A}(\code{"Adj"}),
#' (unnormalized) Laplacian matrix \eqn{L=D-A}(\code{"Lap"}),
#' signless Laplacian matrix \eqn{|L|=D+A}(\code{"SLap"}), or
#' normalized Laplacian matrix \eqn{\tilde{L}=D^{-1/2}LD^{-1/2}}.
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param type type of target structure. One of \code{"Adj","Lap","SLap","NLap"} as defined above.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' \item{spectra}{an \eqn{(N\times M-1)} matrix where each row is top-\eqn{M-1} vibrational spectra.}
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
#' output = nd.dsd(A, out.dist=FALSE)
#' image(output$D, main="two group case")
#'
#' @references
#' \insertRef{wilson_study_2008}{NetworkDistance}
#'
#' @rdname nd_dsd
#' @export
nd.dsd <- function(A, out.dist=TRUE, type=c("Adj","Lap","SLap","NLap")){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.csd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N     = length(listA)    # number of networks
  M     = nrow(listA[[1]]) # number of nodes
  # 3. type : one of 4 matrices
  if (missing(type)){
    type = "Adj"
  } else {
    type = match.arg(type)
  }
  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. setup
  mat_eigs = array(0,c(N,(M-1)))
  mat_dist = array(0,c(N,N))

  #   2. eigenvalue computation
  for (i in 1:N){
    tgt = listA[[i]]
    if (type=="Adj"){
      X = as.matrix(tgt)
    } else if (type=="Lap"){
      X = as.matrix(laplacian_unnormalized(tgt))
    } else if (type=="NLap"){
      X = as.matrix(laplacian_normalized(tgt))
    } else if (type=="SLap"){
      X = as.matrix(laplacian_signless(tgt))
    }
    mat_eigs[i,] = as.vector(RSpectra::eigs(X,(M-1))$values)
  }

  #   3. pairwise computation
  for (i in 1:(N-1)){
    spect1 = mat_eigs[i,]
    for (j in (i+1):N){
      spect2   = mat_eigs[j,]
      solution = sqrt(sum((spect1-spect2)^2))

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
