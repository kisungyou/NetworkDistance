#' 20 adjacency matrices from Erdős–Rényi models
#'
#' Simulated list of 20 adjacency matrices of 28 nodes. First 10 are from Erdős–Rényi model with \eqn{p=0.8}, and
#' the latter 10 are generated using \eqn{p=0.2}. Each element in the list is of size \eqn{(28\times 28)}, symmetric,
#' having values in \eqn{0} or \eqn{1}, and every diagonal element is set as \eqn{0} in accordance with no self-loop assumption.
#'
#' @usage
#' data(graph20)
#'
#' @format
#' A \code{list} of 20 adjacency matrices of size \eqn{(28\times 28)}.
#'
#' @details
#' Below is the code used to generate \emph{graph20}:
#' \preformatted{
#' require(stats)
#' graph20 = list()
#' for (i in 1:10){ # type-1 adjacency matrices
#'   rbin   = rbinom(784,1,0.8)
#'   mat    = matrix(rbin, nrow=28)
#'   matout = mat*t(mat)
#'   diag(matout) = 0
#'   graph20[[i]]=matout
#' }
#' for (i in 11:20){ # type-2 adjacency matrices
#'   rbin   = rbinom(784,1,0.2)
#'   mat    = matrix(rbin, nrow=28)
#'   matout = mat*t(mat)
#'   diag(matout) = 0
#'   graph20[[i]]=matout
#' }
#' }
"graph20"

# graph100 = list()
# for (i in 1:50){ # type-1 adjacency matrices
#  rbin   = rbinom(10000,1,0.8)
#  mat    = matrix(rbin, nrow=100)
#  matout = mat*t(mat)
#  diag(matout) = 0
#  graph100[[i]]=matout
# }
# for (i in 51:100){ # type-2 adjacency matrices
#  rbin   = rbinom(10000,1,0.2)
#  mat    = matrix(rbin, nrow=100)
#  matout = mat*t(mat)
#  diag(matout) = 0
#  graph100[[i]]=matout
# }
