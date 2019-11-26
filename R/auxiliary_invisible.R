# Invisible Auxiliaries ---------------------------------------------------
# 01. graph_transform
#   - support for 'igraph' and 'network' + checking
#   - stop if there is any NA or Inf values conditional on flag
# 02. list_transform
#   - transform the arbitrary data list using 'graph_transform'
#   - check argument : square matrix
# 03. laplacian_unnormalized
# 04. laplacian_normalized
# 05. laplacian_signless
# 06. list_Adj2LapEigs : given a list of adjacency matrices, compute stacked eigenvectors and eigenvalues
#                        mainly used for "gdd" (graph diffusion distance) and laplacian flows.
# 07. list_Adj2LapEigsOrder1 : L1
# 08. list_Adj2SPD           : make laplacian strictly positive semidefinite (global regularization)




# 01. graph_transform -----------------------------------------------------
# sparse Matrix arguments are checked via inherits(obj, "Matrix")
#' @keywords internal
#' @noRd
graph_transform <- function(obj,NIs="allowed"){
  # 1. package:: igraph
  if (inherits(obj, "igraph")){
    output = as_adjacency_matrix(obj)
  # 2. package:: network
  } else if (inherits(obj, "network")){
    output = Matrix(as.matrix.network(obj, matrix.type = "adjacency"), sparse=TRUE)
  # 3. simple matrix
  } else {
    output = Matrix(obj, sparse=TRUE)
  }
  if (NIs!="allowed"){
    if ((any(is.na(output)))||(any(is.infinite(output)))){
      stop("* NetworkDistance : inputs of NA, Inf, or -Inf are not allowed.")
    }
  }
  return(output)
}


# 02. list_transform ------------------------------------------------------
#' @keywords internal
#' @noRd
list_transform <- function(A, NIflag="allowed"){
  # 1. size checker
  n = nrow(A[[1]])
  if (ncol(A[[1]])!=n){
    stop("* NetworkDistance : an input list should contain all square matrices.")
  }
  # 2. transform
  listA = list()
  for (i in 1:length(A)){
    tgt = graph_transform(A[[i]], NIs=NIflag)
    if ((nrow(tgt)!=n)||(ncol(tgt)!=n)){
      stop(paste("* NetworkDistance : ",i,"-st/rd/th matrix in the list has non-matching size.",sep=""))
    }
    listA[[i]] = tgt
  }
  return(listA)
}



# 03. laplacian_unnormalized ----------------------------------------------
#' @keywords internal
#' @noRd
laplacian_unnormalized <- function(matA){
  matD = as.matrix(diag(rowSums(matA))-matA)
  if ((as.double(RSpectra::eigs(matD, 1, which="SM")$values)) < 0){
    matD = as.matrix(Matrix::nearPD(matD, posd.tol=28*.Machine$double.eps)$mat)
  }
  return(matD)
}
# 04. laplacian_normalized ------------------------------------------------
#' @keywords internal
#' @noRd
laplacian_normalized <- function(matA){
  dd = colSums(matA)
  Dinv2 = diag(1/sqrt(dd))
  Dinv2[which(is.infinite(Dinv2))]=0

  D     = diag(dd)
  output = as.matrix(Dinv2%*%(D-matA)%*%Dinv2)
  if ((as.double(RSpectra::eigs(output, 1, which="SM")$values)) < 0){
    output = as.matrix(Matrix::nearPD(output, posd.tol=28*.Machine$double.eps)$mat)
  }
  return(output)
}
# 05. laplacian_signless --------------------------------------------------
#' @keywords internal
#' @noRd
laplacian_signless <- function(matA){
  matD = diag(rowSums(matA))+matA
  return(matD)
}

# 06. list_Adj2LapEigs ----------------------------------------------------
#' @keywords internal
#' @noRd
list_Adj2LapEigs <- function(listA, normalized=FALSE){
  # parameters
  N = length(listA)
  p = nrow(listA[[1]])

  # compute graph laplacians
  Ls = list()
  if (normalized==FALSE){
    for (i in 1:N){
      tgt = listA[[i]]
      diag(tgt) = 0
      Ls[[i]] = laplacian_unnormalized(tgt)
    }
  } else if (normalized==TRUE){
    for (i in 1:N){
      tgt = listA[[i]]
      diag(tgt) = 0
      Ls[[i]] = laplacian_normalized(tgt)
    }
  } else {
    stop("")
  }

  # now do eigendecomposition for each
  eigvals = array(0,c(p,N))   # stack as columns
  eigvecs = array(0,c(p,p,N)) # stack as slices
  for (i in 1:N){
    tgteig = eigen(Ls[[i]])
    eigvals[,i]  = as.vector(tgteig$values)
    eigvecs[,,i] = as.matrix(tgteig$vectors)
  }

  # return
  output = list()
  output$values = eigvals
  output$vectors = eigvecs
  return(output)
}



# 07. list_Adj2LapEigsOrder1 ----------------------------------------------
#' @keywords internal
#' @noRd
list_Adj2LapEigsOrder1 <- function(listA){
  # parameters
  N = length(listA)

  # compute graph laplacians
  Ls = list()
  for (i in 1:N){
    tgt = listA[[i]]
    diag(tgt) = 0
    Ls[[i]] = aux_lapL1(tgt)
  }
  p = nrow(Ls[[i]])

  # now do eigendecomposition for each
  eigvals = array(0,c(p,N))   # stack as columns
  eigvecs = array(0,c(p,p,N)) # stack as slices
  for (i in 1:N){
    tgteig = eigen(Ls[[i]])
    eigvals[,i]  = as.vector(tgteig$values)
    eigvecs[,,i] = as.matrix(tgteig$vectors)
  }

  # return
  output = list()
  output$values = eigvals
  output$vectors = eigvecs
  return(output)
}


# 08. list_Adj2SPD --------------------------------------------------------
#     make laplacian strictly positive semidefinite (global regularization)
#' @keywords internal
#' @noRd
list_Adj2SPD  <- function(listA, normalized=FALSE, stack3d=TRUE){
  # parameters
  N = length(listA)
  p = nrow(listA[[1]])

  # compute graph laplacians
  Ls = list()
  if (normalized==FALSE){
    for (i in 1:N){
      tgt = listA[[i]]
      diag(tgt) = 0
      Ls[[i]] = laplacian_unnormalized(tgt)
    }
  } else if (normalized==TRUE){
    for (i in 1:N){
      tgt = listA[[i]]
      diag(tgt) = 0
      Ls[[i]] = laplacian_normalized(tgt)
    }
  } else {
    stop("")
  }

  # return
  return(Ls)
}

