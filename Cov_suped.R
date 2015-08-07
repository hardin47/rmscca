## Cov.suped returns a (matdim x matdim) covariance matrix for X or Y. ##


Cov.suped <- function(matdim, cor.level, which.matrix, Btype=1, SDX=5){
# Called by sim.setup() and sim.setup.null()
# Calls no other rmscca functions

library(Matrix)


# Create initial D.XX based on "Btype" and "which.matrix".
# D.XX is a (matdim x matdim) diagonal block matrix. Initially, all non-zero entries equal cor.level.
# For each Btype, the variables within a complete group are correlated (as we would expect them to be)
#  For example, if x1 and x2 are in the same block of B, x1 and x2 will be correlated within the simulated data.

if(Btype==0){
    p.dims <- c(2, 3, matdim - 2 - 3)
    q.dims <- c(3, 2, matdim - 3 - 2)
  if(which.matrix=='X'){
    D.XX <- as.matrix(bdiag(matrix(cor.level,p.dims[1],p.dims[1]), matrix(cor.level,p.dims[2],p.dims[2]), 
		matrix(0,p.dims[3],p.dims[3])))
  }else{
    D.XX <- matrix(0, ncol=matdim, nrow=matdim)
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
	var.XX = c(rep((1/cor.level - 1) * ((p.dims[1]^2 - p.dims[1]) * cor.level + p.dims[1]), q.dims[1]), 
		rep((1/cor.level - 1) * ((p.dims[2]^2 - p.dims[2]) * cor.level + p.dims[2]), q.dims[2]), rep(1, q.dims[3]))	
	}
}

# note:  the value for the variance of Y is: cor.level*pX^2 + cor.level*px + 1
# where px is the number of X variables contributing to Y (here 2 and 3 respectively)


if(Btype==1){
    p.dims <- c(12, 8, 12, matdim - 12 - 8 - 12)
    q.dims <- c(8, 4, 4, matdim - 8 - 4 - 4)
  if(which.matrix=='X'){
    D.XX <- as.matrix(bdiag(matrix(cor.level,p.dims[1],p.dims[1]), matrix(cor.level,p.dims[2],p.dims[2]), 
		matrix(cor.level,p.dims[3],p.dims[3]), matrix(0,p.dims[4], p.dims[4])))
  }else{
    D.XX <- matrix(0,ncol=matdim, nrow=matdim)	
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
	var.XX = c(rep((1/cor.level - 1) * ((p.dims[1]^2 - p.dims[1]) * cor.level + p.dims[1]), q.dims[1]), 
		rep((1/cor.level - 1) * ((p.dims[2]^2 - p.dims[2]) * cor.level + p.dims[2]), q.dims[2]),
		rep((1/cor.level - 1) * ((p.dims[3]^2 - p.dims[3]) * cor.level + p.dims[3]), q.dims[3]),
		rep(1, q.dims[4]))
	}
}


if(Btype==2){ 
    p.dims <- c(10, 5, 20, 50, 15, matdim - 10 - 5 - 20 - 50 - 15)
    q.dims <- c(20, 5, 10, 50, 15, matdim - 20 - 5 - 10 - 50 - 15)
  if(which.matrix=='X'){
    D.XX <- as.matrix(bdiag(matrix(cor.level,p.dims[1],p.dims[1]), matrix(cor.level,p.dims[2],p.dims[2]), 
		matrix(cor.level,p.dims[3],p.dims[3]), matrix(cor.level,p.dims[4],p.dims[4]),
		matrix(cor.level,p.dims[5],p.dims[5]), matrix(0,p.dims[6], p.dims[6])))
  }else{
    D.XX <- matrix(0,ncol=matdim, nrow=matdim)	
  }
if(which.matrix=='X'){
	var.XX = rep(1,matdim)     #(matdim x 1) variance vector
  }else{   
	var.XX = c(rep((1/cor.level - 1) * ((p.dims[1]^2 - p.dims[1]) * cor.level + p.dims[1]), q.dims[1]), 
		rep((1/cor.level - 1) * ((p.dims[2]^2 - p.dims[2]) * cor.level + p.dims[2]), q.dims[2]),
		rep((1/cor.level - 1) * ((p.dims[3]^2 - p.dims[3]) * cor.level + p.dims[3]), q.dims[3]),
		rep((1/cor.level - 1) * ((p.dims[4]^2 - p.dims[4]) * cor.level + p.dims[4]), q.dims[4]),
		rep((1/cor.level - 1) * ((p.dims[5]^2 - p.dims[5]) * cor.level + p.dims[5]), q.dims[5]),
		rep(1, q.dims[6]))
	}
}
  
diag(D.XX) <- var.XX		#Each entry of diagonal of D.XX set equal to 1.

# D.XX now created

  

# Pre- and post- multiply D.XX by a standard deviation vector to get a 
#  covariance matrix instead of a correlation matrix.

#return(cov.XX = t( D.XX *  sd.XX) * sd.XX)
return(cov.XX = D.XX)

}
