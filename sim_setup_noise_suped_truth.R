sim.setup.noise.suped.truth <- function(n.obs, B, contamination,var.cor,truth, Btype=1, SDX=1){

  epsilon.val=.01
  edim=25

  library(rrcov)
  library(MASS)
  library(mixstock)
  
  source("Final_funcs/Hardin_Rand_Cov_Func.R")
  source("Final_funcs/sample.sigma12.function.R")
  source("Final_funcs/bic.function.R")
  source("Final_funcs/scca.function.R")
  source("Final_funcs/scca.multiple.R")
  source("Final_funcs/select.parameters.multiple.R")
  source("Final_funcs/randCov_suped.R")
  
  p <- dim(B)[1]
  q <- dim(B)[2]


  firstzero.row = which.min(apply(B,1,sum))
  firstzero.col = which.min(apply(B,2,sum))
  
if(truth == "all"){
  Sigma.g = randCov.suped(p,epsilon.val = epsilon.val,edim = edim,cor.level = var.cor, which.matrix = 'G', Btype, SDX)
  
  G <- mvrnorm(n.obs,rep(0,p), Sigma.g)

  G[,firstzero.row:p] = matrix(rep(0,n.obs*(p-firstzero.row+1)),nrow=n.obs)

  Y.pheno = G%*%B
  
  return(list(X = G, Y = Y.pheno))
}


if(truth == "eps"){

  Sigma.g = randCov.suped(p,epsilon.val = epsilon.val,edim = edim,cor.level = var.cor, which.matrix = 'G', Btype, SDX)
  G <- mvrnorm(n.obs,rep(0,p), Sigma.g)

  Y.pheno = G%*%B
  Sigma.y <- randCov.suped(matdim = q,epsilon.val = epsilon.val, edim = edim, cor.level = var.cor, which.matrix = 'Y',Btype, SDX)
  Y.pheno[,firstzero.col:q] = mvrnorm(n.obs, rep(0,(q-firstzero.col+1)), Sigma.y[firstzero.col:q, firstzero.col:q])

  return(list(X = G, Y = Y.pheno))
}

if(truth == "bigeps"){
  Sigma.g = randCov.suped(p,epsilon.val = epsilon.val,edim = edim,cor.level = var.cor, which.matrix = 'G', Btype, SDX)
  G <- mvrnorm(n.obs,rep(0,p), Sigma.g)
  G[,firstzero.row:p] = mvrnorm(n.obs, rep(0, (p-firstzero.row+1)), 10*Sigma.g[firstzero.row:p, firstzero.row:p])


  mu <- matrix(rep(0,n.obs*q),nrow=n.obs,ncol=q)
  for (i in 1:n.obs){
    mu[i,] <- G[i,]%*%B
  }
  Sigma.y <- randCov.suped(matdim = q,epsilon.val = epsilon.val, edim = edim, cor.level = var.cor, which.matrix = 'Y',Btype, SDX)


  Y.pheno <-list()
  length(Y.pheno) <-q
  for(i in 1:n.obs){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)
    }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  Y.pheno[,firstzero.col:q] = mvrnorm(n.obs, rep(0,(q-firstzero.col+1)), 10*Sigma.y[firstzero.col:q, firstzero.col:q])

  return(list(X = G, Y = Y.pheno))
}




}
