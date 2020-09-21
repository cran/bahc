#' @importFrom fastcluster hclust
#' @import matrixStats
#' @import utils
#' @importFrom stats as.dist
#' @importFrom stats cor rnorm
#'
dendr=function(R){
  out = fastcluster::hclust(as.dist(R),method="average")
  M=out$merge
  #Genealogy Set
  N=nrow(M)

  dend = as.list(1:(N+1))
  names(dend)=as.character(-(1:(N+1)))

  for(i in 1:(N)){
    dend[[as.character(i)]]=list(unlist(dend[[as.character(M[i,1])]]),unlist(dend[[as.character(M[i,2])]]))
  }

  dend=tail(dend,N)
  return(dend)
}

AvLinkC=function(Dend,R){
  N = nrow(R)
  Rs = matrix(0,N,N)

  for(i in seq_along(Dend)){
    a=Dend[[i]][[1]]
    b=Dend[[i]][[2]]
    Rs[a,b] = mean(R[a,b])
  }

  Rs = Rs+t(Rs)
  diag(Rs)=1

  return(Rs)
}


noise=function(N,T,epsilon=1e-10){
  return(matrix(rnorm(N*T,sd = epsilon),nrow=N))
}

#' Compute the BAHC correlation matrix.
#' @export
#' @param x A matrix: \eqn{x_{i,f}} is feature \eqn{f} of object \eqn{i}
#' @param k The order of filtering. \eqn{k=1} corresponds to BAHC.
#' @param Nboot The number of bootstrap copies
#' @return The BAHC-filtered correlation matrix of \code{x}.
#' @examples
#' r=matrix(rnorm(1000),nrow=20)   # 20 objects, 50 features each
#' Cor_bahc=filterCorrelation(r)
filterCorrelation=function(x,k=1,Nboot=100){
  N=dim(x)[1]
  TT=dim(x)[2]

  rT = seq.int(TT)

  Cbav = matrix(0,N,N)
  mynoise = noise(N,TT)
  for(it in 1:Nboot){
    xboot = x[,sample(rT,replace=TRUE)] + mynoise
    Cboot  = cor(t(xboot))
    Cbav = HigherOrder(Cboot,k)
  }
  Cbav = Cbav/Nboot
  return(Cbav)
}


HigherOrder=function(C,k){
  N = dim(C)[1]
  Cf = matrix(0,N,N)
  for(i in 1:k){
    res = C-Cf
    dend = dendr(1-res)
    res = AvLinkC(dend,res)
    diag(res)=0
    Cf = Cf + res
  }
  diag(Cf)= 1
  return(Cf)
}

no_neg=function(C){
  e = eigen(C)
  v = e$vectors
  l = e$values
  mask = l>0
  v = v[,mask,drop=FALSE]
  return( tcrossprod(v * rep(l[mask], each = nrow(v)), v) )
}



#' Compute the BAHC covariance matrix.
#' @export
#' @param x A matrix: \eqn{x_{i,f}} is feature \eqn{f} of object \eqn{i}
#' @param k The order of filtering. \eqn{k=1} corresponds to BAHC.
#' @param Nboot The number of bootstrap copies
#' @return The BAHC-filtered correlation matrix of \code{x}.
#' @examples
#' r=matrix(rnorm(1000),nrow=20)   # 20 objects, 50 features each
#' sigma=exp(runif(20))
#' rs=t(sigma %*% r) %*% sigma
#' Cov_bahc=filterCovariance(rs)
filterCovariance=function(x,k=1,Nboot=100){
  N=dim(x)[1]
  TT=dim(x)[2]

  rT = seq.int(TT)

  Sbav = matrix(0,N,N)
  mynoise=noise(N,TT)
  for(it in 1:Nboot){
    xboot = x[,sample(rT,replace=TRUE)] + mynoise
    Cboot  = cor(t(xboot))
    Cf = HigherOrder(Cboot,k)
    std = rowSds(xboot)
    Sbav = Sbav+ no_neg( Cf*outer(std,std) )
  }
  Sbav = Sbav/Nboot
  return(Sbav)
}
