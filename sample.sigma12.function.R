## Takes X and Y data  and estimates the covariance matrix. ##

sample.sigma12.function <- function(X, Y, robust = T) {
# Takes X and Y data  and estimates the covariance matrix.
# Called by scca.multiple & select.parameters.multiple.R
# Calls no other rmscca functions

# Jake originally had these two lines in the code, but I don't think we want to hard-code
# the training samples
#  x = x.data[training.sample,]
#  y = y.data[training.sample,]
  
# Subtract column (variable) means
X <- X - colMeans(X)
Y <- Y - colMeans(Y)
  
# M-estimation if robust method desired
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
  
# Method From Parkhomenko
# Sample variance-covariance matrices 
sigma11 <- var(X)
sigma22 <- var(Y)

# Change any 0 diagonal entries to 1
diag(sigma11) = ifelse(diag(sigma11)==0, 1, diag(sigma11))
diag(sigma22) = ifelse(diag(sigma22)==0, 1, diag(sigma22))
  
X <- X %*% diag( 1/sqrt(diag(sigma11)) )	#Multiply X by normalized diagonal vector of sigma11
Y <- Y %*% diag( 1/sqrt(diag(sigma22)) )	#Multiply Y by normalized diagonal vector of sigma22
  
if(robust == T){
  return(cov(X,Y, method = "spearman"))

}else{
  return(cov(X,Y,method = "pearson"))
}

}