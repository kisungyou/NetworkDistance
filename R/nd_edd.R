#' Edge Difference Distance
#'
#' It is of the most simplest form that Edge Difference Distance (EDD)
#' takestwo adjacency matrices and takes Frobenius norm of their differnces.
#'
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' }
#'
#' @examples
#' ## generate two types of adjacency matrices of size (3-by-3)
#' rbin1 = rbinom(9,1,0.8); mat1 = matrix(rbin1,nrow=3)
#' rbin2 = rbinom(9,1,0.1); mat2 = matrix(rbin2,nrow=3)
#'
#' mattype1 = ceiling((mat1+t(mat1))/2); diag(mattype1)=0;
#' mattype2 = ceiling((mat2+t(mat2))/2); diag(mattype2)=0;
#'
#' A = list()
#' for (i in 1:3){A[[i]]=mattype1} # first 3 are type-1
#' for (i in 4:6){A[[i]]=mattype2} # next  3 are type-2
#'
#' ## Compute Distance Matrix and Visualize
#' output = nd.edd(A, out.dist=FALSE)
#' image(output$D, main="two group case")
#'
#' @references
#' \insertRef{hammond_graph_2013}{NetworkDistance}
#'
#' @rdname nd_edd
#' @export
nd.edd <- function(A, out.dist=TRUE){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.gdd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  N = length(listA)
  M = nrow(listA[[1]])
  mat_dist = array(0,c(N,N))

  for (i in 1:(N-1)){
    X1 = matrix(listA[[i]], nrow=M)
    for (j in (i+1):N){
      X2 = matrix(listA[[j]], nrow=M)
      solution = aux_FrobeniusDiff(X1, X2)

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
