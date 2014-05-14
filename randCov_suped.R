randCov.suped <- function(matdim, epsilon.val, edim = 25,cor.level, which.matrix, Btype=1, SDX=5){
  library(Matrix)
  epsilon = epsilon.val  # the off-diagonal entries will range from [-epsilon, epsilon]  (not uniformly)
  eidim = edim  # I can tell you more if you want, but we need to select random vectors from some big dimension
  
  #1.Generate a diagonal matrix with 1s on the diagonal (then subtract epsilon)
 # D.XX = diag(rep(1 - epsilon,matdim))  # a pxp matrix with 1 - epsilon on the diagonal.
  

if(Btype==0){
  if(which.matrix=='G'){
  D.XX <- as.matrix(bdiag(matrix(rep(cor.level,4),2),matrix(rep(cor.level,9),3),
		matrix(rep(0,(matdim-2-3)^2),(matdim-2-3))))
  }else{
    D.XX <- as.matrix(bdiag(matrix(rep(cor.level,9),3),matrix(rep(cor.level,4),2),
		matrix(rep(0,(matdim-3-2)^2),(matdim-3-2))))
  }
  }
if(Btype==1){
  if(which.matrix=='G'){
  D.XX <- as.matrix(bdiag(matrix(rep(cor.level,144),12),matrix(rep(cor.level,64),8), matrix(rep(cor.level,144),12),
		matrix(rep(0,(matdim-8-12-12)^2),(matdim-8-12-12))))
  }else{
    D.XX <- as.matrix(bdiag(matrix(rep(cor.level,64),8),matrix(rep(cor.level,16),4), matrix(rep(cor.level,16),4),

		matrix(rep(0,(matdim-8-4-4)^2),(matdim-8-4-4))))
  }
  }

if(Btype==2){  # needs to be fixed
  if(which.matrix=='G'){
  D.XX <- as.matrix(bdiag(matrix(rep(cor.level,144),12),matrix(rep(cor.level,64),8), matrix(rep(cor.level,144),12),
		matrix(rep(0,(matdim-8-12-12)^2),(matdim-8-12-12))))
  }else{
    D.XX <- as.matrix(bdiag(matrix(rep(cor.level,64),8),matrix(rep(cor.level,16),4), matrix(rep(cor.level,16),4),
		matrix(rep(0,(matdim-8-4-4)^2),(matdim-8-4-4))))
  }
  }

  
  diag(D.XX) <- (rep(1-epsilon,matdim))
    
  #2.  Add random entries between -epsilon and epsilon such that the matrix
  #   remains positive definite (and therefore remains a correlation matrix).
  
  eivect = c()
  for ( i in 1:matdim){
    ei = runif(eidim, -1, 1)
    eivect = cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2)))
  }
  bigE = t(eivect) %*% eivect
  cor.XX = D.XX + bigE
  
  #3 .Pre- and post- multiply the whole thing by a st dev vector to get a 
  #   covariandce matrix instead of a correlation matrix.
  
  sd.XX = rnorm(matdim, SDX, 1)  

# I don’t know how you want to come up with the variances on the diagonal…  Make sure they are all positive.
  
  return(cov.XX = t( cor.XX *  sd.XX) * sd.XX)}
