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
#' ## load example data
#' data(graph20)
#'
#' ## Compute Distance Matrix and Visualize
#' output = nd.edd(graph20, out.dist=FALSE)
#' opar   = par(pty="s")
#' image(output$D[,20:1], main="two group case", axes=FALSE, col=gray(0:32/32))
#' par(opar)
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
