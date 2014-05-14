sim.setup <- function(n.obs, B){
  library(rrcov)
  library(MASS)
  library(mixstock)
  
  source("Current_R_Files/Hardin's_Random_Covariance_Function.R")
  source("My_Functions2/sample.sigma12.function.R")
  source("My_Functions2/bic.function.R")
  source("My_Functions2/scca.function.R")
  source("My_Functions2/scca.multiple.R")
  source("My_Functions2/select.parameters.multiple.R")
  
  p <- dim(B)[1]
  q <- dim(B)[2]
  
  Sigma.g = randCov(p,epsilon.val = .5,edim = 25)
  G.clean <- mvrnorm(n=n.obs, mu = rep(0,p),Sigma = Sigma.g)
  
  mu <- matrix(rep(0,n.obs*q),nrow=n.obs,ncol=q)
  for (i in 1:n.obs){
    mu[i,] <- G.clean[i,]%*%B
  }
  
  #need to find reasonable parameters for this
  Sigma.y <- randCov(matdim = q,epsilon.val = .5, edim = 25)
  
  Y.pheno <-list()
  length(Y.pheno) <-q
  for(i in 1:n.obs){
    Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)
  }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  return(list(X = G.clean, Y = Y.pheno))
}
