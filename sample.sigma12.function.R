sample.sigma12.function <- function(x, y, robust = T)
{

# Jake originally had these two lines in the code, but I don't think we want to hard-code
# the training samples
#  x = x.data[training.sample,]
#  y = y.data[training.sample,]
  
  x <- t(t(x) - apply(x,2,mean))
  y <- t(t(y) - apply(y,2,mean))
  
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

  diag(sigma11) = ifelse(diag(sigma11)==0,1,diag(sigma11))
  diag(sigma22) = ifelse(diag(sigma22)==0,1,diag(sigma22))
  
  x <- x %*% diag( 1/sqrt(diag(sigma11)) )
  y <- y %*% diag( 1/sqrt(diag(sigma22)) )
  
  if(robust == T){
    return(cov(x,y, method = "spearman"))
  }else{
    return(cov(x,y,method = "pearson"))
  }
}
