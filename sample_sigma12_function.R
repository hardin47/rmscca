## Takes X and Y data  and estimates the covariance matrix. ##

sample.sigma12.function <- function(X, Y, robust = TRUE) {
# Takes X and Y data  and estimates the covariance matrix.
# Called by scca.multiple & select.parameters.multiple.R
# Calls no other rmscca functions
  
# Subtract column (variable) means
X <- X - colMeans(X)
Y <- Y - colMeans(Y)
  
  
# Method From Parkhomenko
# Sample variance-covariance matrices 
sigma11 <- var(X)
sigma22 <- var(Y)

# Change any 0 diagonal entries to 1
diag(sigma11) = ifelse(diag(sigma11)==0, 1, diag(sigma11))
diag(sigma22) = ifelse(diag(sigma22)==0, 1, diag(sigma22))
  
X <- X %*% diag( 1/sqrt(diag(sigma11)) )	#Multiply X by normalized diagonal vector of sigma11
Y <- Y %*% diag( 1/sqrt(diag(sigma22)) )	#Multiply Y by normalized diagonal vector of sigma22
  
if(robust == TRUE){
  return(cov(X,Y, method = "spearman"))

}else{
  return(cov(X,Y,method = "pearson"))
}

}