#' @keywords internal
#' @noRd
aux_lapL1 <- function(A){
  E<-sum(A)/2 # sum of adjacent matrix gives twice the number of edges
  if (E==1){
    return(array(0,c(1,1)))
  } else {
    #by the hand shake lemma



    #Define entries sum function ES, which sums the upper triangle
    #entries of the adjacent matrix to (ij) position, which enumerate the edges
    ES<-function(p,q){
      result<-0
      if (q<=p)
      {result<-NA}
      else if(p==1)
        result<-sum(A[p,1:q])
      else
      {
        for(i in 1:(p-1))
        {result<-result+sum(A[i,(i+1):ncol(A)])}
        result<-result+sum(A[p,(p+1):q])
      }
      return(result)
    }

    #compute A_l the lower adjacent matrix


    Al<-matrix(0,nrow=E,ncol=E)

    B<-list()                                ###########creat a list, the nth element is a list [ij] containing the nodes

    ############## of n-th edge of the graph G
    for(i in 1:(ncol(A)-1))
    {for(j in (i+1):ncol(A))
    {if (A[i,j]==1)
    {n<-ES(i,j)
    B[[n]]<-c(i,j)}
    }
    }

    #define the entries of Al, the lower adjacent matrix

    for(i in 1:(E-1))
    {for(j in (i+1):E)
    {if((B[[i]][1]==B[[j]][1])|(B[[i]][2]==B[[j]][1])|(B[[i]][2]==B[[j]][2]))
    {Al[i,j]<-1
    Al[j,i]=Al[i,j]
    }
    }

    }
    #This finishes the lower adjacent matrix




    #start programming degree matrix



    A2<-A%*%A                   ################## A2 is the square of the adjacent matrix
    D1<-matrix(0,nrow=E,ncol=E)
    for(r in 1:(ncol(A)-1))
    {for(s in (r+1):ncol(A))                        #define the degree the the nth Edge $ij$
    { if (A[r,s]==1)                            #using the $ij$ entry of the square
    {m<-ES(r,s)                                 #of the adjacent matrix A
    D1[m,m]<-A2[r,s]
    }
    }
    }

    #start programming upper adjacent matrix Au

    Au<-matrix(0,nrow=E,ncol=E)
    for(i in 1:(E-1))
    {for(j in (i+1):E)
    {if(((B[[i]][1]==B[[j]][1])&(A[B[[i]][2],B[[j]][2]]==1))|((B[[i]][2]==B[[j]][1])&(A[B[[i]][1],B[[j]][2]]==1))|((B[[i]][2]==B[[j]][2])&(A[B[[i]][1],B[[j]][1]]==1)))
    {Au[i,j]<-1
    Au[j,i]<-1
    }
    }
    }

    #Define the L1 Laplacian operator L1:= D1-Au+Al+2*I

    L1<-D1-Au+Al+2*diag(E)
    return(L1)
  }
}
