#######################################
#-------------SCCA For Multiple Lambdas----------------#
####################################

scca.multiple <- function(X,Y, lambda.u = rep(0,dim(X)[2]), lambda.v = rep(0,dim(Y)[2]), rob = T, bic = F){
  p <- dim(X)[2] #number variables X
  q <- dim(Y)[2] #number variables Y
  n.sample<-dim(X)[1] #number samples in X and Y
  num.pairs <- min(p,q)
  
  lambda.u <- numeric(num.pairs)
  lambda.v <- numeric(num.pairs)
  
  coefs.u.list <- list()
  coefs.v.list <- list()
  nobic.coef.u = c()
  nobic.coef.v = c()
  coef.u = c()
  coef.v = c()
  u.initial <- list()
  v.initial <- list()
  bic.cor <- c()
  nobic.cor <- c()
  
  length(coefs.u.list) <- num.pairs
  length(coefs.v.list) <- num.pairs
  length(u.initial) <- num.pairs
  length(u.initial) <- num.pairs
  
  #First coef pairs are calculated
  #to avoid SNPs with variance of 0, we divide by the total sample sigma rather than the testing sample sigma
  sigma11.sample <- var(X)
  k <- list()
  length(k) <- num.pairs
  k[[1]] <- sample.sigma12.function(X,Y, robust = rob)#, sigma11.sample) #parameter 'rob' dictates whether or not robust method 
  
  #Initial vectors, as chosen by Parkhomenko et al.
  u.initial.first <- as.numeric(k[[1]] %*% rep(1,q)/q)
  u.initial[[1]] <- u.initial.first/sqrt(as.numeric(t(u.initial.first)%*%u.initial.first))
  v.initial.first <- as.numeric(t(k[[1]]) %*% rep(1,p)/p)
  v.initial[[1]] <- v.initial.first /sqrt(as.numeric(t(v.initial.first)%*%v.initial.first))
  
  #Use CV to select parameters if either is set to 0
  if(lambda.u[1] == 0 | lambda.v[1] == 0){
    parameters <- select.parameters.multiple(X,Y, bound = 2, count = 1, rob.p = rob)
    lambda.u[1] <- parameters$lambda.u[1]
    lambda.v[1] <- parameters$lambda.v[1]
  }
  
  
  uv <- list()
  length(uv) <- num.pairs
  
  #Finding the first sparse vectors
  uv[[1]] <- scca.function(k[[1]], u.initial[[1]], v.initial[[1]], lambda.u[1], lambda.v[1])
  
  
  coefs.u.list[[1]]<-uv[[1]]$u.new
  coefs.v.list[[1]]<-uv[[1]]$v.new
  
  var.x = var(X)
  var.y = var(Y)
  diag(var.x) = ifelse(diag(var.x) ==0, 1, diag(var.x))
  diag(var.y) = ifelse(diag(var.y) ==0, 1, diag(var.y))


  
  d <- numeric(num.pairs)
  d[1] <- as.numeric(t(coefs.u.list[[1]])%*%k[[1]]%*%coefs.v.list[[1]]) #coefficient for decomposition
  if(rob==T){
    nobic.cor[1] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[1]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[1]], method = "spearman")
  }else{
    nobic.cor[1] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[1]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[1]], method = "pearson")
  }
  
  #Calculating the remaining sparse vectors
  for(i in 2:num.pairs){
    k[[i]]<- k[[i-1]] - d[i-1]*coefs.u.list[[i-1]]%*%t(coefs.v.list[[i-1]]) #residual matrix
    
    u.initial.first <- as.numeric(k[[i]] %*% rep(1,q)/q)
    u.initial[[i]] <- u.initial.first/sqrt(as.numeric(t(u.initial.first)%*%u.initial.first))
    v.initial.first <- as.numeric(t(k[[i]]) %*% rep(1,p)/p)
    v.initial[[i]] <- v.initial.first /sqrt(as.numeric(t(v.initial.first)%*%v.initial.first))
    
    parameters <- select.parameters.multiple(X,Y, bound = 2,n.cv = 5, d.vec = d,u.list =  coefs.u.list, v.list =  coefs.v.list, count = i, rob.p = rob)
    lambda.u[i] <- parameters$lambda.u[1]
    lambda.v[i] <- parameters$lambda.v[1]
    
    
    uv[[i]] <- scca.function(k[[i]], u.initial[[i]], v.initial[[i]], lambda.u[i], lambda.v[i])
    
    
    coefs.u.list[[i]]<-uv[[i]]$u.new
    coefs.v.list[[i]]<-uv[[i]]$v.new
    
    
    d[i] <- as.numeric(t(coefs.u.list[[i]])%*%k[[i]]%*%coefs.v.list[[i]]) #coefficient for decomposition
  
    #new
    if(rob==T){
      nobic.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[i]], method = "spearman")
    }else{
      nobic.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[i]], method = "pearson")
    }
  }
  
  #new
  nobic.coef.u <- data.frame(matrix(unlist(coefs.u.list), nrow=length(coefs.u.list[[1]]), byrow=F))
  nobic.coef.v <- data.frame(matrix(unlist(coefs.v.list), nrow=length(coefs.v.list[[1]]), byrow=F))
  
  
  #BIC - what if we did this after decomposition?
  if(bic ==T){
    for(i in 1:num.pairs){
      if(abs(d[i])>0){
        coefficients <- bic.function(uv[[i]],X,Y,p,q,n.sample,robust=rob)
        coefs.u.list[[i]]<-coefficients$coefs.u
        coefs.v.list[[i]]<-coefficients$coefs.v
      } 
      if(rob==T){
        bic.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[i]], method = "spearman")
      }else{
        bic.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%coefs.v.list[[i]], method = "pearson")
      }
    }
    coef.u <- data.frame(matrix(unlist(coefs.u.list), nrow=length(coefs.u.list[[1]]), byrow=F))
    coef.v <- data.frame(matrix(unlist(coefs.v.list), nrow=length(coefs.v.list[[1]]), byrow=F))
    
  }
  
  #turning lists into data frames 
 
  return(list(nobic.coef.u = nobic.coef.u, nobic.coef.v = nobic.coef.v, bic.coef.u = coef.u, bic.coef.v =coef.v, bic.cor = bic.cor, nobic.cor = nobic.cor, lambda.u = lambda.u, lambda.v = lambda.v))
}
