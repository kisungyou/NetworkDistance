#' Network Flow Distance
#'
#' @param A a list of length \eqn{N} containing adjacency matrices.
#' @param order the order of Laplacian; currently only 0 and 1 are supported.
#' @param out.dist a logical; \code{TRUE} for computed distance matrix as a \code{dist} object.
#' @param vect a vector of parameters \eqn{t} whose values will be used.
#'
#' @return a named list containing \describe{
#' \item{D}{an \eqn{(N\times N)} matrix or \code{dist} object containing pairwise distance measures.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(23)
#' Total<-20
#' N1<-Total/2
#' P1<-0.75
#' P2<-0.6
#' P12=0.04
#' Iteration<-2
#' CAP<-4
#'
#' bb<-list()         ## edges to remove
#' bb[[1]]<-c(1,1)
#' bb[[2]]<-c(4,19)
#' bb[[3]]<-c(12,17)
#' bb[[4]]<-c(13,18)
#' bb[[5]]<-c(1,3)
#' bb[[6]]<-c(15,8)
#' bb[[7]]<-c(2,6)
#'
#'
#' A<-matrix(0,nrow=Total,ncol=Total)                       ######### define adjacent matrix
#' for(i in (1:(N1-1)))
#' {for(j in ((i+1):N1))
#' {A[i,j]<-rbinom(1,1,P1)
#' A[j,i]<-A[i,j]}
#' }
#'
#' for(i in ((N1+1):(Total-1)))
#' {for(j in ((i+1):Total))
#' { A[i,j]<-rbinom(1,1,P2)
#' A[j,i]<-A[i,j]
#' }
#' }
#' for(i in (1:N1))
#' {for(j in (N1+1):Total)
#' {A[i,j]<-rbinom(1,1,P12)
#' A[j,i]<-A[i,j]
#' }
#' }
#'
#' listA = list()
#' for (i in 1:7){
#'    tgtA = A
#'    idm  = bb[[i]][1]
#'    idn  = bb[[i]][2]
#'
#'    tgtA[idm,idn] = 0
#'    tgtA[idn,idm] = 0
#'    listA[[i]] = tgtA
#' }
#'
#' # compute two diffusion-based distances and visualize
#' out1 = nd.gdd(listA, out.dist=FALSE)$D
#' out2 = testdec(listA, out.dist=FALSE)$D
#' par(mfrow=c(1,2))
#' image(pracma::flipud(out1),col=gray((0:32)/32), main="Hammond Pairwise Distance",axes=FALSE)
#' image(pracma::flipud(out2),col=gray((0:32)/32), main="Dianbin Pairwise Distance",axes=FALSE)
#' }
#'
#' @export
nd.nfd <- function(A, order=0, out.dist=TRUE, vect=seq(from=0,to=10,length.out=1000)){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.nfd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")
  # 3. vect
  if ((!is.vector(vect))||(any(vect<0))||(any(is.na(vect)))||(any(is.infinite(vect)))){
    stop("* nd.nfd : input 'vect' should be a vector of nonnegative real numbers.")
  }
  vect = sort(vect)

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  # N = length(listA)
  # mat_dist = array(0,c(N,N))
  #
  # for (i in 1:(N-1)){
  #   L1 = as.matrix(gdd_laplacian(listA[[i]]))
  #   for (j in (i+1):N){
  #     L2    = as.matrix(gdd_laplacian(listA[[j]]))
  #
  #     L12dist = lfdistance(L1,L2,0.1)
  #     mat_dist[i,j] = L12dist
  #     mat_dist[j,i] = L12dist
  #   }
  # }

  N = length(listA)
  if (order==0){
    Lprocess = list_Adj2LapEigs(listA)
  } else if (order==1){
    if (length(unique(unlist(lapply(listA, sum))))!=1){
      stop("* nd.nfd : for the order 1 case, all networks must have same number of edges.")
    }
    Lprocess = list_Adj2LapEigsOrder1(listA)
  } else {
    stop("* nd.nfd : orders other than k=0,1 are not supported.")
  }
  Lvecs    = Lprocess$vectors
  Lvals    = Lprocess$values
  mat_dist = array(0,c(N,N))

  for (i in 1:(N-1)){
    L1 = as.matrix(Lvecs[,,i])
    D1 = as.vector(Lvals[,i])
    for (j in (i+1):N){
      L2 = as.matrix(Lvecs[,,j])
      D2 = as.vector(Lvals[,j])

      distvalue = lfdistance_new(L1,D1,L2,D2,vect)
      mat_dist[i,j] = distvalue
      mat_dist[j,i] = distvalue
    }
  }


  #-------------------------------------------------------
  ## Return output
  if (out.dist){
    mat_dist = as.dist(mat_dist)
  }
  output = list()
  output$D = mat_dist
  return(output)
}
