## Simulates the data to test ##

sim.setup.noise.suped <- function(n.obs, B, contamination, var.cor, noisetype, Btype=1, SDX=1){

# Called by big_sim_cutoff.R
# Calls Cov.suped
# source("Final_funcs/Cov_suped.R")
# library(MASS)
  
p = dim(B)[1]
q = dim(B)[2]
  
if(noisetype == "clean"){
  Sigma.x = Cov.suped(p,cor.level = var.cor, which.matrix = 'X', Btype, SDX=1)	#Create X covariance matrix
  Sigma.y = Cov.suped(q, cor.level = var.cor, which.matrix = 'Y',Btype, SDX=1)	#Create Y covariance matrix
  
  # Produces n.obs samples from multivariate normal distibrution with mean 0 and covariance matrix Sigma.x
  X = mvrnorm(n.obs,rep(0,p), Sigma.x)		
  mu = X%*%B


  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in 1:n.obs){
  	Y.pheno[i,] = mvrnorm(1, mu[i,], Sigma.y)
#	Y.pheno[i,] = mvrnorm(1, rep(0,q), Sigma.y)
  }

  return(list(X = X, Y = Y.pheno))
}



if(noisetype=="t"){

  tdf = 2
  Sigma.x = Cov.suped(p,cor.level = var.cor, which.matrix = 'X', Btype, SDX)	#Create X covariance matrix
  Sigma.y = Cov.suped(q, cor.level = var.cor, which.matrix = 'Y',Btype, SDX)	#Create Y covariance matrix
  
  # Produces n.obs samples from multivariate t-distibrution with mean 0 and covariance matrix Sigma.x
  X = mvrnorm(n.obs,rep(0,p), Sigma.x)/sqrt(rchisq(n.obs,tdf)/tdf)	
  mu = X%*%B

  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in 1:n.obs){
  	Y.pheno[i,] = mvrnorm(1, mu[i,], Sigma.y)/sqrt(rchisq(1,tdf)/tdf)
  }

  return(list(X = X, Y = Y.pheno))
}



if(noisetype=="sym"){

  Sigma.x = Cov.suped(p, cor.level = var.cor, which.matrix = 'X', Btype, SDX)
  Sigma.y = Cov.suped(matdim = q, cor.level = var.cor, which.matrix = 'Y', Btype, SDX)

  X = mvrnorm(n.obs, rep(0,p), Sigma.x)
  mu = X %*%B

  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in 1:n.obs){
  	Y.pheno[i,] = mvrnorm(1, mu[i,], Sigma.y)
  }

  mins.X = min(apply(X,1,min))
  maxes.X = max(apply(X,1,max))
  mins.Y = min(apply(Y.pheno,1,min))
  maxes.Y = max(apply(Y.pheno,1,max))

  for(i in 1:n.obs){
  	if(runif(1) < contamination){
		    X[i,] = runif(p, 3*mins.X, 3*maxes.X)
		    Y.pheno[i,] = runif(q, 3*mins.Y, 3*maxes.Y)
		    }
	}


  return(list(X = X, Y = Y.pheno))

}


if(noisetype=="asym"){

  Sigma.x = Cov.suped(p, cor.level = var.cor, which.matrix = 'X', Btype, SDX)
  Sigma.y = Cov.suped(matdim = q, cor.level = var.cor, which.matrix = 'Y', Btype, SDX)

  num.noise <- n.obs*contamination
  num.clean <- n.obs*(1-contamination)

  X1 <- jitter(matrix(mean(diag(Sigma.x)),num.noise,p))  # noise
  X2 <- mvrnorm(num.clean, rep(0,p), Sigma.x)		# clean
  X.full <- rbind(X1,X2)				# full

  ourorder <- sample(nrow(X.full))
  X.full <- X.full[ourorder,]				# scramble the order of noise
 
  mu = X.full%*%B

  outliers <- ourorder[ourorder<=num.noise]
  nonoutliers <- ourorder[ourorder>num.noise]

  Y.pheno = matrix(ncol=q, nrow=n.obs)
  for(i in outliers){
      Y.pheno[i,] <- jitter(rep(sum(diag(Sigma.y)), q))
     }
  for(i in nonoutliers){
      Y.pheno[i,] <- mvrnorm(1,mu[i,], Sigma.y)
  }
  return(list(X = X.full, Y = Y.pheno))

}



}
