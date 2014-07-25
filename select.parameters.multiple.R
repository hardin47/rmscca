## Figures out which lambda values are best to use for scca on X and Y ##

select.parameters.multiple <- function(X, Y, bound = 2, n.cv = 5, d.vec, u.list, v.list,  count, rob.p = T){
# Called by scca.multiple
# Calls sample.sigma12.function() & scca.function



lambda.v.seq <- seq(0, bound, by=0.1)  # Possible values of sparseness parameters for data Y. Lower bounds should be 0, upper bound can be increased to 2.
lambda.u.seq <- seq(0, bound, by=0.1)  # Possible values of sparseness parameters for data X. Lower bounds should be 0, upper bound can be increased to 2.

n.lambdas.v <-  length(lambda.v.seq) 	# Number of possible sparseness parameters for data Y.
n.lambdas.u <-  length(lambda.u.seq)	# Number of possible sparseness parameters for data X.

lambda.v.matrix <- matrix(rep(lambda.v.seq, n.lambdas.u), nrow=n.lambdas.u, byrow=T)
lambda.u.matrix <- matrix(rep(lambda.u.seq, n.lambdas.v), nrow=n.lambdas.u, byrow=F)

p <- dim(X)[2]			  	# Number of X variables
q <- dim(Y)[2]			 	# Number of Y variables
n.obs <- dim(X)[1]			# Number of observations


ones.q <- rep(1, q)/q			# Vector length q, entries = 1/q
ones.p <- rep(1, p)/p			# Vector length p, entries = 1/p


###_______________________________________________________________________________________
###Analysis
###_______________________________________________________________________________________

n.cv.sample <- trunc(n.obs/n.cv)	# Number of observations in training sample (trunc will round down if n.obs/n.cv is not an integer
whole.sample <- seq(1, n.obs)		# Vector length n.obs, all entries = 1. (Markers for determining testing/training samples?)


predict.corr.scca <- matrix(0, nrow=n.lambdas.u, ncol=n.lambdas.v)		# Average test sample correlation for each combination of
								        	#  sparseness parameters will be entered into predict.corr.scca
   

#_______Cross-validation to select optimal combination of sparseness parameters_______#

for (i.cv in 1:n.cv) 	# n.cv is number of cross-validations
{
  testing.sample <- whole.sample[((i.cv-1)*n.cv.sample+1):(i.cv*n.cv.sample)]	# Determining testing sample (subset of whole.sample) elements, of size n.cv.sample
  training.sample <- whole.sample[!whole.sample%in%testing.sample]		# Determining training sample elements, all of whole.sample except testing sample


  # Estimate covariance matrix of training simple
  k <- sample.sigma12.function(X[training.sample, ], Y[training.sample, ], robust=rob.p)


  # Calculating the residual matrix  (note that for the (count)th pair, we remove (count-1) pieces of information)
  if(count>1){
    for(l in 2:count){
      k <- k - d.vec[l-1]*u.list[[l-1]]%*%t(v.list[[l-1]])
    }
  }



  # Get starting values for singular vectors as column and row means from matrix K
  u.initial <- k %*% ones.q						# Step (2,a,i) from Coleman: u.initial =  alpha = K*Beta (better way to say this...)
  u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))	# Normalize u.initial
  v.initial <- t(k) %*% ones.p						# Step (2,b,i) from Coleman: v.initial = beta = K*alpha   (better way to say this...)
  v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial)) 	# Normalize v.initial


  #_______________Data for Predicted correlation (testing sample)_______________#

  x.predict <- X[testing.sample, ]			# X data of testing sample
  y.predict <- Y[testing.sample, ]			# Y data of testing sample

  # Standardize data
  x.predict <- t(t(x.predict)-apply(x.predict,2,mean))
  y.predict <- t(t(y.predict)-apply(y.predict,2,mean))

  sigma11.predict <- var(x.predict)
  sigma22.predict <- var(y.predict)


  # Make sure the diagonal variance entries aren't zero, if so, change to 1.
  diag(sigma11.predict) = ifelse(diag(sigma11.predict)==0,1,diag(sigma11.predict))
  diag(sigma22.predict) = ifelse(diag(sigma22.predict)==0,1,diag(sigma22.predict))

  x.predict <- x.predict %*% diag( 1/sqrt(diag(sigma11.predict)) )
  y.predict <- y.predict %*% diag( 1/sqrt(diag(sigma22.predict)) )


  #____________Loops for sparseness parameter combinations____________#
  for(j.lambda.v in 1:n.lambdas.v)
  { 
    for(j.lambda.u in 1:n.lambdas.u) 		# All combinations of lambdas tested on u.inital and v.initial
    {
      lambda.v <- lambda.v.seq[j.lambda.v]	# sparseness parameter for Y
      lambda.u <- lambda.u.seq[j.lambda.u]	# sparseness parameter for X

      uv <- scca.function(k, u.initial, v.initial, lambda.u, lambda.v)

      vj <- uv$v.new
      uj <- uv$u.new

      #print(c(j.lambda.v,j.lambda.u))

      # Calculate predicted correlation for SCCA
      if(rob.p == T){
        predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + (cor(x.predict%*%uj, y.predict%*%vj, method = "spearman"))
      }else{
        predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + (cor(x.predict%*%uj, y.predict%*%vj, method = "pearson"))
      }

      if(is.na(predict.corr.scca[j.lambda.u, j.lambda.v])) {  # close if
	predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] <-  NA
        break 
      }
  
    }	# close loop on lambda.u
  }	# close loop on lambda.v
    
}  # close cross-validation loop
  
  
#______________Identify optimal sparseness parameter combination______________#
predict.corr.scca[is.na(predict.corr.scca)] <- 0
predict.corr.scca <- predict.corr.scca/n.cv

best.lambda.v <- lambda.v.matrix[which(predict.corr.scca==max(predict.corr.scca,na.rm=T))]
best.lambda.u <- lambda.u.matrix[which(predict.corr.scca==max(predict.corr.scca,na.rm=T))]

return(list(best.cor = best.predict.corr.scca, lambda.v = best.lambda.v, lambda.u = best.lambda.u))

}