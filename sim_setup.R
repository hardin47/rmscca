## Simulates the data to test ##

sim.setup <- function(n.obs, B, contamination, var.cor, noisetype, Btype=1, SDX=1){

# Called by fullSimScript.R
# Calls Cov.suped
# source("Final_funcs/Cov_suped.R")
# library(mvnfast)
  
p = dim(B)[1]
q = dim(B)[2]
  
if(noisetype == "clean"){
  Sigma.x = Cov.suped(p,cor.level = var.cor, which.matrix = 'X', Btype, SDX=1)	#Create X covariance matrix
  Sigma.y = Cov.suped(q, cor.level = var.cor, which.matrix = 'Y',Btype, SDX=1)	#Create Y covariance matrix
  
  # Produces n.obs samples from multivariate normal distibrution with mean 0 and covariance matrix Sigma.x
  X = rmvn(n.obs,rep(0,p), Sigma.x)		
  mu = X%*%B


  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in 1:n.obs){
  	Y.pheno[i,] = rmvn(1, mu[i,], Sigma.y)
  }



  return(list(X = X, Y = Y.pheno))
}



if(noisetype=="t"){

  tdf = 2
  Sigma.x = Cov.suped(p,cor.level = var.cor, which.matrix = 'X', Btype, SDX)	#Create X covariance matrix
  Sigma.y = Cov.suped(q, cor.level = var.cor, which.matrix = 'Y',Btype, SDX)	#Create Y covariance matrix
  
  # Produces n.obs samples from multivariate t-distibrution with mean 0 and covariance matrix Sigma.x
  X = rmvn(n.obs,rep(0,p), Sigma.x)/sqrt(rchisq(n.obs,tdf)/tdf)	
  mu = X%*%B

  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in 1:n.obs){
  	Y.pheno[i,] = rmvn(1, mu[i,], Sigma.y)/sqrt(rchisq(1,tdf)/tdf)
  }

  return(list(X = X, Y = Y.pheno))
}

}
