#' tester
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
#'}
#'
#' @export
testdec <- function(A){
  #-------------------------------------------------------
  ## PREPROCESSING
  # 1. list of length larger than 1
  if ((!is.list(A))||(length(A)<=1)){
    stop("* nd.csd : input 'A' should be a list of length larger than 1.")
  }
  # 2. transform the data while checking
  listA = list_transform(A, NIflag="not")

  #-------------------------------------------------------
  ## MAIN COMPUTATION
  N = length(listA)
  mat_dist = array(0,c(N,N))

  for (i in 1:(N-1)){
    L1 = as.matrix(gdd_laplacian(listA[[i]]))
    for (j in (i+1):N){
      L2    = as.matrix(gdd_laplacian(listA[[j]]))

      L12dist = lfdistance(L1,L2,0.1)
      mat_dist[i,j] = L12dist
      mat_dist[j,i] = L12dist
    }
  }

  return(mat_dist)
}
