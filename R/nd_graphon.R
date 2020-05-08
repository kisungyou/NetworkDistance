#' Graphon Estimates Distance
#'
#' Graphon is a symmetric measurable function
#' \deqn{W:[0,1]^2\rightarrow[0,1]}
#' that is considered to be a generating model for an observed network. \code{nd.graphon} computes
#' distances between networks based on the estimated graphons of each network. Estimation methods
#' are taken from \pkg{graphon} package. For more details, see the function links below.
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param method type of graphon estimation methods to be used.
#' @param ... extra parameters to be passed onto graphon estimation functions. See also \code{\link[graphon]{est.completion}},
#' \code{\link[graphon]{est.LG}}, \code{\link[graphon]{est.nbdsmooth}}, \code{\link[graphon]{est.SBA}}, and \code{\link[graphon]{est.USVT}} for details.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' }
#'
#' @examples
#' ## load example data
#' data(graph20)
#'
#' ## compute USVT-based distance
#' output <- nd.graphon(graph20, out.dist=FALSE, method="usvt")
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(output$D[,20:1], main="USVT", col=gray(0:32/32), axes=FALSE)
#' par(opar)
#'
#' @references
#' \insertRef{mukherjee_clustering_2017}{NetworkDistance}
#'
#' @rdname nd_graphon
#' @export
nd.graphon <- function(A, out.dist=TRUE, method=c("completion","LG","nbdsmooth","SBA","USVT"), ...){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.graphon : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N     = length(listA)    # number of graphs
  p     = nrow(listA[[1]]) # number of nodes
  # 3. method
  allmethods = tolower(c("completion","LG","nbdsmooth","SBA","USVT"))
  if (missing(method)){
    mymethod = "usvt"
  } else {
    mymethod = match.arg(tolower(method), allmethods)
  }
  # 4. return type
  myreturn = as.logical(out.dist)

  #-------------------------------------------------------
  # Main Computation
  # 1. compute graphon
  listG = list()
  for (i in 1:N){
    tgtA = as.matrix(listA[[i]])  # select the target
    if (all(mymethod=="completion")){
      listG[[i]] = graphon::est.completion(A=tgtA, ...)
    } else if (all(mymethod=="lg")){
      listG[[i]] = graphon::est.LG(A=tgtA, ...)$P
    } else if (all(mymethod=="nbdsmooth")){
      listG[[i]] = graphon::est.nbdsmooth(A=tgtA)$P
    } else if (all(mymethod=="sba")){
      listG[[i]] = graphon::est.SBA(A=tgtA, ...)$P
    } else if (all(mymethod=="usvt")){
      listG[[i]] = graphon::est.USVT(A=tgtA, ...)$P
    }
  }
  # 2. parwise distance
  mat_dist = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      mat_dist[i,j] = base::norm(listG[[i]] - listG[[j]], type="F")
      mat_dist[j,i] = mat_dist[i,j]
    }
  }

  #-------------------------------------------------------
  ## RETURN RESULTS
  if (myreturn){
    mat_dist = as.dist(mat_dist)
  }

  result = list()
  result$D= mat_dist
  return(result)
}


# ## load example data
# data(graph20)
#
# ## compute USVT-based distance
# out1 <- nd.graphon(graph20, out.dist=FALSE, method="completion")
# out2 <- nd.graphon(graph20, out.dist=FALSE, method="LG")
# out3 <- nd.graphon(graph20, out.dist=FALSE, method="nbdsmooth")
# out4 <- nd.graphon(graph20, out.dist=FALSE, method="SBA")
# out5 <- nd.graphon(graph20, out.dist=FALSE, method="USVT")
#
# ## visualize
# opar = par(no.readonly=TRUE)
# par(mfrow=c(2,3), pty="s")
# image(out1$D[,20:1], col=gray(0:32/32), axes=FALSE, main="completion")
# image(out2$D[,20:1], col=gray(0:32/32), axes=FALSE, main="LG")
# image(out3$D[,20:1], col=gray(0:32/32), axes=FALSE, main="nbdsmooth")
# image(out4$D[,20:1], col=gray(0:32/32), axes=FALSE, main="SBA")
# image(out5$D[,20:1], col=gray(0:32/32), axes=FALSE, main="USVT")
# par(opar)
