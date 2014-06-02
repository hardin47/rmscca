Cov.suped <- function(matdim, cor.level, which.matrix, Btype=1, SDX=5){

  library(Matrix)
    

# for each Btype, the variables within a complete group are correlated (as we would expect them to be)

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

  
  diag(D.XX) <- rep(1,matdim)
    
  #   Pre- and post- multiply the whole thing by a st dev vector to get a 
  #   covariandce matrix instead of a correlation matrix.
  
  sd.XX = rnorm(matdim, SDX, 1)  
  cor.XX = D.XX

# I don’t know how you want to come up with the variances on the diagonal…  Make sure they are all positive.
  
  return(cov.XX = t( cor.XX *  sd.XX) * sd.XX)}
