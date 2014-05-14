randCov <- function(matdim = 100, epsilon.val = .9, edim = 25, SDX = 5){
  p=matdim  # ??
  epsilon = epsilon.val  # the off-diagonal entries will range from [-epsilon, epsilon]  (not uniformly)
  eidim = edim  # I can tell you more if you want, but we need to select random vectors from some big dimension
  
  #1.Generate a diagonal matrix with 1s on the diagonal (then subtract epsilon)
  D.XX = diag(rep(1 - epsilon,p))  # a pxp matrix with 1 - epsilon on the diagonal.
  
  #2.  Add random entries between -epsilon and epsilon such that the matrix
  #   remains positive definite (and therefore remains a correlation matrix).
  
  eivect = c()
  for ( i in 1:p){
    ei = runif(eidim, -1, 1)
    eivect = cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2)))
  }
  bigE = t(eivect) %*% eivect
  cor.XX = D.XX + bigE
  
  #3 .Pre- and post- multiply the whole thing by a st dev vector to get a 

  sd.XX=rnorm(p, SDX, 1)

  #   covariance matrix instead of a correlation matrix.
  return(cov.XX = t( cor.XX *  sd.XX) * sd.XX)}