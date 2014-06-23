## Simulates the data to test ##

sim.setup.noise.suped <- function(n.obs, B, contamination, var.cor, noisetype, Btype=1, SDX=1){
# Called by big_sim_cutoff.R
# Calls Cov.suped
  
source("Final_funcs/Cov_suped.R")
  
p <- dim(B)[1]
q <- dim(B)[2]
  
if(noisetype == "clean"){
  Sigma.g = Cov.suped(p,cor.level = var.cor, which.matrix = 'G', Btype, SDX)	#Create X covariance matrix
  
  G <- mvrnorm(n.obs,rep(0,p), Sigma.g)		#Produces n.obs samples from multivariate distibrution
						# with means 0 and covariance matrix Sigma.g
  mu <- matrix(0,nrow=n.obs,ncol=q)
  mu <- G%*%B
  
  # Need to find reasonable parameters for this.
  Sigma.y <- Cov.suped(q, cor.level = var.cor, which.matrix = 'Y',Btype, SDX)	#Create Y covariance matrix

  Y.pheno <-list()
  length(Y.pheno) <-q

  for(i in 1:n.obs){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)
    }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  return(list(X = G, Y = Y.pheno))
}



if(noisetype=="t"){

  tdf = 1
  Sigma.g = Cov.suped(p, cor.level = var.cor, which.matrix = 'G', Btype, SDX)

  G <- mvrnorm(n.obs,rep(0,p), Sigma.g)/sqrt(rchisq(n.obs,tdf)/tdf)
  mu <- matrix(rep(0,n.obs*q),nrow=n.obs,ncol=q)
  for (i in 1:n.obs){
    mu[i,] <- G[i,]%*%B
  }

  #need to find reasonable parameters for this
  Sigma.y <- Cov.suped(matdim = q, cor.level = var.cor, which.matrix = 'Y', Btype, SDX)

  Y.pheno <-list()
  length(Y.pheno) <-q
  for(i in 1:n.obs){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)/sqrt(rchisq(q,tdf)/tdf)
    }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  return(list(X = G, Y = Y.pheno))

}


if(noisetype=="sym"){

  Sigma.g = Cov.suped(p, cor.level = var.cor, which.matrix = 'G', Btype, SDX)

  #need to find reasonable parameters for this
  Sigma.y <- Cov.suped(matdim = q, cor.level = var.cor, which.matrix = 'Y', Btype, SDX)

  G.clean = mvrnorm(n.obs, rep(0,p), Sigma.g)

  Y.pheno = G.clean%*%B

#  mu = G.clean %*%B
#  Y.pheno = matrix(ncol=q, nrow=n.obs)
#  for(i in 1:n.obs){
#  	Y.pheno[i,] = mvrnorm(1, mu[i,], Sigma.y)
#	}

  mins.G = min(apply(G.clean,1,min))
  maxes.G = max(apply(G.clean,1,max))
  mins.Y = min(apply(Y.pheno,1,min))
  maxes.Y = max(apply(Y.pheno,1,max))

  for(i in 1:n.obs){
  	if(runif(1) < contamination){
		    G.clean[i,] = runif(p, 3*mins.G, 3*maxes.G)
		    Y.pheno[i,] = runif(q, 3*mins.Y, 3*maxes.Y)
		    }
	}


  return(list(X = G.clean, Y = Y.pheno))

}


if(noisetype=="asym"){

  Sigma.g = Cov.suped(p,  cor.level = var.cor, which.matrix = 'G', Btype, SDX)

  num.noise <- n.obs*contamination
  num.clean <- n.obs*(1-contamination)

  G1 <- jitter(matrix(rep(mean(diag(Sigma.g)),num.noise*p), ncol=p))
  G2 <- mvrnorm(num.clean, rep(0,p), Sigma.g)
  G.clean <- rbind(G1,G2)

  ourorder <- sample(nrow(G.clean))

  G.clean <- G.clean[ourorder,]

  G.clean <- matrix(unlist(G.clean),nrow = n.obs, byrow =T)

  mu <- matrix(rep(0,n.obs*q),nrow=n.obs,ncol=q)
  for (i in 1:n.obs){
    mu[i,] <- G.clean[i,]%*%B
  }

  #need to find reasonable parameters for this
  Sigma.y <- Cov.suped(matdim = q,  cor.level = var.cor, which.matrix = 'Y', Btype, SDX)

  outliers <- ourorder[ourorder<=num.noise]
  nonoutliers <- ourorder[ourorder>num.noise]

  Y.pheno <-list()
  length(Y.pheno) <-q
  for(i in outliers){
      Y.pheno[[i]] <- rep(sum(diag(Sigma.y)), q)
     }
  for(i in nonoutliers){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)
  }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  return(list(X = G.clean, Y = Y.pheno))

}



}
