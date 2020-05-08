#' Log Moments Distance
#'
#' For a graph with an adjacency matrix \eqn{A}, \emph{graph moment} is defined as
#' \deqn{\rho_m (A) = tr(A/n)^m}
#' where \eqn{n} is the number of vertices and \eqn{m} is an order of the moment. \code{nd.moments} computes
#' pairwise distances based on log of graph moments from \eqn{m=1} to \eqn{m=k}.
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param k the integer order of moments. If \eqn{k} is too large, it may incur numerical overflow.
#' @param metric type of distance measures for log-moment features. See \code{\link[stats]{dist}} for more details.
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
#' ## compute distance based on different k's.
#' out3 <- nd.moments(graph20, k=3, out.dist=FALSE)
#' out5 <- nd.moments(graph20, k=5, out.dist=FALSE)
#' out7 <- nd.moments(graph20, k=7, out.dist=FALSE)
#' out9 <- nd.moments(graph20, k=9, out.dist=FALSE)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(out3$D[,20:1], col=gray(0:32/32), axes=FALSE, main="k=3")
#' image(out5$D[,20:1], col=gray(0:32/32), axes=FALSE, main="k=5")
#' image(out7$D[,20:1], col=gray(0:32/32), axes=FALSE, main="k=7")
#' image(out9$D[,20:1], col=gray(0:32/32), axes=FALSE, main="k=9")
#' par(opar)
#'
#' @references
#' \insertRef{mukherjee_clustering_2017}{NetworkDistance}
#'
#' @rdname nd_moments
#' @export
nd.moments <- function(A, k=3, metric=c("euclidean","maximum","manhattan","canberra","binary","minkowski"), out.dist=TRUE){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.moments : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N     = length(listA)    # number of graphs
  # 3. k and metric
  myk = round(k)
  if (missing(metric)){
    mymetric = "euclidean"
  } else {
    mymetric = match.arg(metric)
  }
  # 4. return type
  myreturn = as.logical(out.dist)

  #-------------------------------------------------------
  # MAIN COMPUTATION
  # 1. compute features
  features = array(0,c(N,myk))
  for (n in 1:N){
    features[n,] = (moment_single2(as.matrix(listA[[n]]), myk))
  }
  # 2. compute distance
  dmat = stats::dist(features, method=mymetric)

  #-------------------------------------------------------
  ## RETURN RESULTS
  if (!myreturn){
    dmat = as.matrix(dmat)
  }

  result = list()
  result$D= dmat
  return(result)
}


# single ------------------------------------------------------------------
# requires Amat to be
# - single1 follows Somendu's code.
# - single2 with the definition. I go with definition.
#' @keywords internal
#' @noRd
moment_single1 <- function(Amat, k){
  n = nrow(Amat)
  temp = Amat
  mult = Amat

  output = rep(0,k)
  for (i in 1:k){
    output[i] = (base::sum(base::diag(temp)))/(exp(1+(i/2.0)*log(n)))
    if (i < k){
      temp = temp%*%mult
    }
  }
  return(log(output))
}
#' @keywords internal
#' @noRd
moment_single2 <- function(Amat, k){
  n = nrow(Amat)
  temp = Amat/n
  mult = Amat/n

  output = rep(0,k)
  for (i in 1:k){
    output[i] = base::sum(base::diag(temp))
    temp = temp%*%mult
  }
  return(log(output))
}
