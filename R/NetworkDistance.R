#' Distance Measures for Networks
#'
#' Network has gathered much attention from many disciplines, as many of real data
#' can be well represented in the relational form. The concept of distance - or, metric - between
#' two networks is the starting point for inference on population of networks. \pkg{NetworkDistance} package
#' provides a not-so-comprehensive collection of distance measures for measuring dissimilarity between two network objects.
#' Data should be supplied as \emph{adjacency} matrices, where we support three formats of data representation;
#' \code{matrix} object in \pkg{R} base, \code{network} class from \pkg{network} package, and \code{igraph} class from
#' \pkg{igraph} package.
#'
#' @docType package
#' @name NetworkDistance
#' @import Rdpack
#' @import Matrix
#' @import RSpectra
#' @import CovTools
#' @importFrom utils packageVersion
#' @importFrom pracma flipud
#' @importFrom igraph as_adjacency_matrix graph_from_adjacency_matrix degree closeness betweenness
#' @importFrom network as.matrix.network
#' @importFrom stats as.dist integrate rbinom
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom parallel detectCores stopCluster makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom Rcpp evalCpp
#' @useDynLib NetworkDistance
NULL

