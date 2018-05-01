#' HIM Distance
#'
#' Hamming-Ipsen-Mikhailov (HIM) combines the local Hamming edit distance and the global
#' Ipsen-Mikhailov distance to merge information at each scale. For Ipsen-Mikhailove distance,
#' it is provided as \code{nd.csd} in our package for consistency. Given a parameter \eqn{\xi} (\code{xi}),
#' it is defined as
#' \deqn{HIM_{\xi}(A,B)=\sqrt{H^2(A,B)+\xi\cdot IM^2(A,B)}/\sqrt{1+\xi}}
#' where \eqn{H} and \eqn{IM} stand for Hamming and I-M distance, respectively.
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param xi a parameter to control balance between two distances.
#' @param ntest the number of searching over \code{\link{nd.csd}} parameter.
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
#' ## compute distance and visualize
#' output = nd.him(A, out.dist=FALSE)
#' image(output$D, main="two group case")
#'
#' @references
#' \insertRef{jurman_him_2015}{NetworkDistance}
#'
#' @seealso \code{\link{nd.hamming}}, \code{\link{nd.csd}}
#' @rdname nd_him
#' @export
nd.him <- function(A, out.dist=TRUE, xi=1.0, ntest=10){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.him : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  N = length(listA)
  M = nrow(listA[[1]])
  # 3. xi : parameter
  xi = as.double(xi)
  if ((length(as.vector(xi))!=1)||(xi<0)||(is.na(xi))||(is.infinite(xi))){
    stop("* nd.him : parameter 'xi' should be a nonnegative real number.")
  }
  # 4. ntest : parameter
  ntest = as.integer(ntest)
  if ((length(as.vector(ntest))!=1)||(ntest<2)||(is.na(ntest))||(is.infinite(ntest))){
    stop("* nd.him : parameter 'ntest' should be a positive integer larger than 1.")
  }

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. compute Hamming Distance
  D_hamming <- nd.hamming(listA, out.dist=FALSE)$D
  #   2. IM-type
  #   2-1. we need to find parameter,  the optimal one
  mat_empty = matrix(rep(0,M*M),nrow=M)
  mat_full  = matrix(rep(1,M*M),nrow=M); diag(mat_full)=0;
  A_IMfind  = list()
  A_IMfind[[1]] = mat_empty; A_IMfind[[2]] = mat_full;

  im_grid  = 10^seq(from=-2,to=1,length.out=ntest)
  im_score = rep(0,ntest)
  for (i in 1:ntest){
    im_scoremat <- nd.csd(A_IMfind, out.dist=FALSE, bandwidth = im_grid[i])$D
    im_score    <- im_scoremat[1,2]
  }
  im_index = which(abs(im_score-1)==(min(abs(im_score - 1))))
  im_optbd = im_grid[im_index]

  #   2-2. let's compute the distance matrix
  D_im <- nd.csd(listA, out.dist=FALSE, bandwidth = im_optbd)$D

  #   3. now we have to matrices, go compute it.
  solution = (sqrt((D_hamming^2)+ xi*(D_im^2)))/sqrt(1+xi)

  #-------------------------------------------------------
  ## RETURN RESULTS
  if (out.dist){
    solution = as.dist(solution)
  }

  result = list()
  result$D= solution
  return(result)
}
