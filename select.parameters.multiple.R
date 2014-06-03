select.parameters.multiple <- function(x.data, y.data, bound = 2, n.cv = 5, d.vec, u.list, v.list,  count, rob.p = T){
  
  lambda.v.seq <- seq(0, bound, by=0.02)  # Possible values of sparseness parameters for data Y. Lower bouns should be 0, upper bound can be increased to 2.
  lambda.u.seq <- seq(0, bound, by=0.02)  # Possible values of sparseness parameters for data X. Lower bouns should be 0, upper bound can be increased to 2.
  
  n.lambdas.u <-  length(lambda.u.seq)
  n.lambdas.v <-  length(lambda.v.seq)
  
  lambda.v.matrix <- matrix(rep(lambda.v.seq, n.lambdas.u), nrow=n.lambdas.u, byrow=T)
  lambda.u.matrix <- matrix(rep(lambda.u.seq, n.lambdas.v), nrow=n.lambdas.u, byrow=F)
  
  p <- dim(x.data)[2]
  q <- dim(y.data)[2]
  n.sample <- dim(x.data)[1]
  
  ones.p <- rep(1, p)/p
  ones.q <- rep(1, q)/q
  
  
  ### _______________________________________________________________________________________
  ### Analysis
  ### _______________________________________________________________________________________
  
  n.cv.sample <- trunc(n.sample/n.cv)
  whole.sample <- seq(1, n.sample)
  
  predict.corr.scca <- matrix(0, nrow=n.lambdas.u, ncol=n.lambdas.v)  # This matrix will contain average test sample correlation for each combination of sparseness parameters
  
  
  #_______Cross-validation to select optimal combination of sparseness parameters____________
  for (i.cv in 1:n.cv)
  {
    testing.sample <- whole.sample[((i.cv-1)*n.cv.sample+1):(i.cv*n.cv.sample)]
    training.sample <- whole.sample[!whole.sample%in%testing.sample]
    
    k <- sample.sigma12.function(x.data[training.sample, ], y.data[training.sample, ], robust=rob.p)
    
    #Calculating the residual matrix  (note that for the (count)th pair, we remove (count-1) pieces of information)
    if(count>1){
      for(l in 2:count){
        k <- k - d.vec[l-1]*u.list[[l-1]]%*%t(v.list[[l-1]])
      }
    }
    
    # Get starting values for singular vectors
    # as column and row means from matrix K
    u.initial <- k %*% ones.q
    u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
    v.initial <- t(k) %*% ones.p
    v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))
    
    # _______________Data for Predicted correlation (testing sample)_________________
    
    x.predict <- x.data[testing.sample, ]
    y.predict <- y.data[testing.sample, ]
    
    # Standardize data
    x.predict <- t(t(x.predict)-apply(x.predict,2,mean))
    y.predict <- t(t(y.predict)-apply(y.predict,2,mean))
    
    sigma11.predict <- var(x.predict)
    sigma22.predict <- var(y.predict)


    # making sure the diagonal variance entries aren't zero
    diag(sigma11.predict) = ifelse(diag(sigma11.predict)==0,1,diag(sigma11.predict))
    diag(sigma22.predict) = ifelse(diag(sigma22.predict)==0,1,diag(sigma22.predict))
    
    x.predict <- x.predict %*% diag( 1/sqrt(diag(sigma11.predict)) )
    y.predict <- y.predict %*% diag( 1/sqrt(diag(sigma22.predict)) )
    
    
    # ____________Loops for sparseness parameter combinations__________
    for(j.lambda.v in 1:n.lambdas.v)
    {
      
      flag.na <- 0
      
      for(j.lambda.u in 1:n.lambdas.u)
      {
        lambda.v <- lambda.v.seq[j.lambda.v]	# sparseness parameter for Y
        lambda.u <- lambda.u.seq[j.lambda.u]	# sparseness parameter for X
        
        if(flag.na==0)
        {
          uv <- scca.function(k, u.initial, v.initial, lambda.u, lambda.v)
          
          vj <- uv$v.new
          uj <- uv$u.new

	#print(c(j.lambda.v,j.lambda.u))
          
          # Calculate predicted correlation for SCCA
          if(rob.p == T){
            predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + abs(cor(x.predict%*%uj, y.predict%*%vj, method = "spearman"))
          }else{
            predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + abs(cor(x.predict%*%uj, y.predict%*%vj, method = "pearson"))
          }
          
          if(is.na(predict.corr.scca[j.lambda.u, j.lambda.v])) flag.na <- 1
        }	# close if
        
        if(flag.na==1)
        {
          predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] <- predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] + NA
          break
        }
        
      }	# close loop on lambda.u
    }	# close loop on lambda.v
    
  }	# close cross-validation loop
  
  
  # ______________Identify optimal sparseness parameter combination___________		
  
  predict.corr.scca[is.na(predict.corr.scca)] <- 0
  predict.corr.scca <- predict.corr.scca/n.cv
  
  best.predict.corr.scca <- max(abs(predict.corr.scca), na.rm=T)
  best.lambda.v <- lambda.v.matrix[predict.corr.scca==best.predict.corr.scca]
  best.lambda.u <- lambda.u.matrix[predict.corr.scca==best.predict.corr.scca]
  
  return(list(best.cor = best.predict.corr.scca, lambda.v = best.lambda.v, 
              lambda.u = best.lambda.u))
}
