sim.setup.noise <- function(n.obs, B, contamination){
  library(rrcov)
  library(MASS)
  library(mixstock)
  
  source("Current_R_Files/Hardin's_Random_Covariance_Function.R")
  source("My_Functions/sample.sigma12.function.R")
  source("My_Functions/bic.function.R")
  source("My_Functions/scca.function.R")
  source("My_Functions/scca.multiple.R")
  source("My_Functions/select.parameters.multiple.R")
  
  p <- dim(B)[1]
  q <- dim(B)[2]
  
  Sigma.g = randCov(p,epsilon.val = .9,edim = 25)
  
#   outliers <- which(runif(n.obs,0,100)<5)
#   nonoutliers <- setdiff(1:n.obs,outliers)

#   G.clean <- list()
#   length(G.clean) <- p
#   for(i in outliers){
#       G.clean[[i]] <- mvrnorm(1,rep(0,p),9*Sigma.g)
#     }
#   for(i in nonoutliers){
#       G.clean[[i]] <- mvrnorm(1,rep(0,p),Sigma.g)
#   }
  num.noise <- n.obs*contamination
  num.clean <- n.obs*(1-contamination)
  
  G1 <- mvrnorm(num.noise,rep(0,p), 9*Sigma.g)
  G2 <- mvrnorm(num.clean, rep(0,p), Sigma.g)
  G.clean <- rbind(G1,G2)
  
  order <- sample(nrow(G.clean))
  
  G.clean <- G.clean[order,]
  
  G.clean <- matrix(unlist(G.clean),nrow = n.obs, byrow =T)
  
  mu <- matrix(rep(0,n.obs*q),nrow=n.obs,ncol=q)
  for (i in 1:n.obs){
    mu[i,] <- G.clean[i,]%*%B
  }
  
  #need to find reasonable parameters for this
  Sigma.y <- randCov(matdim = q,epsilon.val = .9, edim = 25)
  
  outliers <- order[order<=num.noise]
  nonoutliers <- order[order>num.noise]
  
  Y.pheno <-list()
  length(Y.pheno) <-q
  for(i in outliers){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], 9*Sigma.y)
    }
  for(i in nonoutliers){
      Y.pheno[[i]] <- mvrnorm(1,mu[i,], Sigma.y)
  }
  Y.pheno <- matrix(unlist(Y.pheno),nrow = n.obs, byrow =T)
  return(list(X = G.clean, Y = Y.pheno))
}
