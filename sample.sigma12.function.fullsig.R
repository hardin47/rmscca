sample.sigma12.function.fullsig <- function(x, y, x.sub, y.sub, robust = T)
{
  x <- t(t(x) - apply(x,2,mean))
  y <- t(t(y) - apply(y,2,mean))

  x.sub = t(t(x.sub) - apply(x.sub,2,mean))
  y.sub = t(t(y.sub) - apply(y.sub,2,mean))
  
  #M-estimation if robust method desired
  #     Big = cbind(x,y)
  #     r.cov = CovMest(Big)@cov
  #     dimx = dim(x)[2]
  #     dimy = dim(y)[2]
  #     dims = dimy+dimx
  #     Sigma11 <- diag(1/sqrt(diag(r.cov[1:dimx,1:dimx])))
  #     
  #     Sigma22 <- diag(1/sqrt(diag(r.cov[(dimx+1):dims, (dimx+1):dims])))
  #     Sigma12 <- r.cov[1:dimx,(dimx+1):dims]
  #     
  #     return(Sigma11%*%Sigma12%*%Sigma22)
  
  #Method From Parkhomenko
  # Sample variance-covariance matrices 
  sigma11 <- var(x)
  sigma22 <- var(y)

  
  x.sub <- x.sub %*% diag( 1/sqrt(diag(sigma11)) )
  y.sub <- y.sub %*% diag( 1/sqrt(diag(sigma22)) )
  
  if(robust == T){
    return(cov(x.sub,y.sub, method = "spearman"))
  }else{
    return(cov(x.sub,y.sub,method = "pearson"))
  }
}
