#' \eqn{L_2} Distance of Continuous Spectral Densities
#'
#' The method employs spectral density of eigenvalues from
#' Laplacian in that for each, we have corresponding
#' spectral density \eqn{\rho(w)} as a sum of
#' narrow Lorentz distributions with \code{bandwidth} parameter.
#' Since it involves integration of a function over the
#' non-compact domain, it may blow up to infinity and the code
#' automatically aborts the process.
#'
#' @examples
#' \donttest{
#' ## load example data
#' data(graph20)
#'
#' ## compute distance matrix
#' output = nd.csd(graph20, out.dist=FALSE, bandwidth=1.0)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(output$D[,20:1], main="two group case", axes=FALSE, col=gray(0:32/32))
#' par(opar)
#' }
#'
#' @param A a list of length \eqn{N} containing \eqn{(M\times M)} adjacency matrices.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param bandwidth common bandwidth of positive real number.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' \item{spectra}{an \eqn{(N\times M-1)} matrix where each row is top-\eqn{M-1} vibrational spectra.}
#' }
#'
#' @references
#' \insertRef{ipsen_evolutionary_2002}{NetworkDistance}
#'
#' @rdname nd_csd
#' @export
nd.csd <- function(A, out.dist=TRUE, bandwidth=1.0){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.csd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  # 3. bandwidth : positive bandwidth
  bandwidth = as.double(bandwidth)
  if ((length(as.vector(bandwidth))!=1)||(bandwidth<=0)||(is.na(bandwidth))||(is.infinite(bandwidth))){
    stop("* nd.csd : parameter 'bandwidth' should be a positive real number.")
  }

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. setup
  N = length(listA)
  M = nrow(listA[[1]])
  mat_seig = list()
  mat_dist = array(0,c(N,N))
  vec_Ks   = rep(0,N)

  #   2. eigenvalue computations at once + transform into Vibrational Spectra
  for (i in 1:N){
    L             = as.matrix(laplacian_unnormalized(listA[[i]]))
    mat_seig[[i]] = sqrt(as.vector(RSpectra::eigs(L,(M-1))$values))
    vec_Ks[i]     = tryCatch(compute_normc(mat_seig[[i]], bandwidth),
                             error=function(e){stop("* csd : non-finite normalizing constant. please use other methods.")})
  }

  #   3. pairwise computation
  for (i in 1:(N-1)){
    spect1 = mat_seig[[i]]
    K1     = vec_Ks[i]
    for (j in (i+1):N){
      spect2 = mat_seig[[j]]
      K2     = vec_Ks[2]

      integrand = function(w, k1, k2, sp1, sp2, bd){
        N = length(sp1); output = 0;
        for (i in 1:N){
          output = output + (k1*(w/(((w-sp1[i])^2)+(bd^2)))) - (k2*(w/(((w-sp2[i])^2)+(bd^2))))
        }
        return((output^2))
      }
      sol = stats::integrate(integrand, lower=0, upper=Inf, subdivisions=1234,
                             k1=K1, k2=K2, sp1=spect1, sp2=spect2, bd=bandwidth)$value
      if ((is.na(sol))||(is.infinite(sol))){
        solution = NA
      } else {
        solution = sqrt(sol)
      }
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
  result$spectra = mat_seig
  return(result)
}


#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
compute_normc <- function(sp, bd){
  Kinv = stats::integrate(function(w,bd,sp){output=0;N=length(sp);
  for (i in 1:N){
    output = output + (w/(((w-sp[i])^2)+(bd^2)))
  }
  return(output)}, lower = 0, upper = Inf, bd=bd, sp=sp, subdivisions=1234)
  return(1/Kinv$value)
}
