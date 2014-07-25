## Cov.suped returns a (matdim x matdim) covariance matrix for X or Y. ##


Cov.suped <- function(matdim, cor.level, which.matrix, Btype=1, SDX=5){
# Called by sim.setup.noise.suped()
# Calls no other rmscca functions

library(Matrix)


# Create initial D.XX based on "Btype" and "which.matrix".
# D.XX is a (matdim x matdim) diagonal block matrix. Initially, all non-zero entries equal cor.level.
# For each Btype, the variables within a complete group are correlated (as we would expect them to be)
#  For example, if x1 and x2 are in the same block of B, x1 and x2 will be correlated within the simulated data.

if(Btype==0){
  if(which.matrix=='X'){
    D.XX <- as.matrix(bdiag(matrix(cor.level,2,2), matrix(cor.level,3,3), matrix(0,matdim-2-3,matdim-2-3)))
  }else{
    D.XX <- as.matrix(bdiag(matrix(cor.level,3,3), matrix(cor.level,2,2), matrix(0,matdim-3-2,matdim-3-2)))
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
	var.XX = c(rep(4,2), rep(7,3), rep(1,matdim-3-2))	
	}
}

# note:  the value for the variance of Y is: var.cor*pX^2 + var.cor*px + 1
# where px is the number of X variables contributing to Y (here 2 and 3 respectively)


if(Btype==1){
  if(which.matrix=='X'){
    D.XX <- as.matrix(bdiag(matrix(cor.level,12,12), matrix(cor.level,8,8), matrix(cor.level,12,12), 	
	  	 	 matrix(0,matdim-8-12-12,matdim-8-12-12)))
  }else{
    D.XX <- as.matrix(bdiag(matrix(cor.level,8,8), matrix(cor.level,4,4), matrix(cor.level,4,4),
			matrix(0,matdim-8-4-4,matdim-8-4-4)))	
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
	var.XX = c(rep(79,8), rep(37,4), rep(79,4), rep(1,matdim-8-4-4))	
	}
}


if(Btype==2){ 
  if(which.matrix=='X'){
  D.XX <- as.matrix(bdiag(matrix(cor.level,10,10),matrix(cor.level,5,5), matrix(cor.level,20,20),
		matrix(cor.level,50,50), matrix(cor.level,15,15),
		matrix(0,matdim-10-5-20-50-15,matdim-10-5-20-50-15)))
  }else{
    D.XX <- as.matrix(bdiag(matrix(cor.level,20,20),matrix(cor.level,5,5), matrix(cor.level,10,10),
		matrix(cor.level,50,50), matrix(cor.level,15,15),
		matrix(0,matdim-20-5-10-50-15,matdim-20-5-10-50-15)))
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
#	var.XX = c(rep(56,20), rep(16,5), rep(211,10), rep(1276,50), rep(121,15), rep(1,matdim-20-5-10-50-15))
	var.XX = rep(1,matdim)
	}
}
  
diag(D.XX) <- var.XX		#Each entry of diagonal of D.XX set equal to 1.

# D.XX now created

  

# Pre- and post- multiply D.XX by a standard deviation vector to get a 
#  covariance matrix instead of a correlation matrix.

#return(cov.XX = t( D.XX *  sd.XX) * sd.XX)
return(cov.XX = D.XX)

}
