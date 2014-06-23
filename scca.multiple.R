#####################################################
#-------------SCCA For Multiple Lambdas-------------#
#####################################################

scca.multiple <- function(X, Y, rob = T) {
# Called by big_sim_cutoff.R
# Calls sample.sigma12.function(), select.parameters.multiple() & scca.function()


p <- dim(X)[2] 	# Number of variables X
q <- dim(Y)[2]	# Number of variables Y
minpq <- min(p,q)	

# Initialize sparseness parameters vectors and left/right singular vector lists
lambda.u <- numeric(minpq)
lambda.v <- numeric(minpq)
u.initial <- list()
v.initial <- list()

sp.coefs.u.list <- list()
sp.coefs.v.list <- list()
sp.coef.u <- c()
sp.coef.v <- c()
sp.cor <- c()

length(sp.coefs.u.list) <- minpq
length(sp.coefs.v.list) <- minpq
length(u.initial) <- minpq
length(u.initial) <- minpq

####  SPARSE CCA  ####
##  We use CV to find the right lambda.u and lambda.v values.
##  Because we are finding multiple canonical pairs, we start with the first pair, then
##   subtract out the signal from each pair before finding the subsequent pair  

  
# First coeficient pairs are calculated
# To avoid SNPs with variance of 0, we divide by the total sample sigma rather than the testing sample sigma
sigma11.sample <- var(X)
k <- list()
length(k) <- minpq
k[[1]] <- sample.sigma12.function(X,Y, robust = rob) # Parameter 'rob' dictates whether or not robust method 


# Initial singular vectors: Column means of k for v or row means of K for u, standardized to unit length. (Parkhomenko)
u.initial.first <- as.numeric(k[[1]] %*% rep(1,q)/q)
u.initial[[1]] <- u.initial.first/sqrt(as.numeric(t(u.initial.first)%*%u.initial.first))
v.initial.first <- as.numeric(t(k[[1]]) %*% rep(1,p)/p)
v.initial[[1]] <- v.initial.first /sqrt(as.numeric(t(v.initial.first)%*%v.initial.first))


# Use CV to select parameters if either is set to 0
# The select.parameters.multiple function finds the best lambda values for the given X and Y
if(lambda.u[1] == 0 | lambda.v[1] == 0){
  parameters <- select.parameters.multiple(X,Y, bound = 1 count = 1, rob.p = rob)
  lambda.u[1] <- parameters$lambda.u[1]
  lambda.v[1] <- parameters$lambda.v[1]
}


uv <- list()
length(uv) <- minpq


# Finding the first sparse vectors of coefficients (u.new and v.new)
uv[[1]] <- scca.function(k[[1]], u.initial[[1]], v.initial[[1]], lambda.u[1], lambda.v[1])


# Save the first pair of canonical coefficients
sp.coefs.u.list[[1]]<-uv[[1]]$u.new
sp.coefs.v.list[[1]]<-uv[[1]]$v.new


var.x = var(X)
var.y = var(Y)
diag(var.x) = ifelse(diag(var.x) ==0, 1, diag(var.x))
diag(var.y) = ifelse(diag(var.y) ==0, 1, diag(var.y))


if(rob==T){
  sp.cor[1] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%sp.coefs.u.list[[1]],Y%*%diag(sqrt(1/diag(var.y)))%*%sp.coefs.v.list[[1]], method = "spearman")
}else{
  sp.cor[1] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%sp.coefs.u.list[[1]],Y%*%diag(sqrt(1/diag(var.y)))%*%sp.coefs.v.list[[1]], method = "pearson")
}


# Break down K so that we can find subsequent canonical pairs
# Coefficients for decomposition
d <- numeric(minpq)
d[1] <- as.numeric(t(sp.coefs.u.list[[1]])%*%k[[1]]%*%sp.coefs.v.list[[1]])


# Calculating the remaining sparse vectors
for(i in 2:minpq){

  # Break down K so that we can find subsequent canonical pairs (residual matrix)
  k[[i]]<- k[[i-1]] - d[i-1]*sp.coefs.u.list[[i-1]]%*%t(sp.coefs.v.list[[i-1]])

  u.initial.first <- as.numeric(k[[i]] %*% rep(1,q)/q)
  u.initial[[i]] <- u.initial.first/sqrt(as.numeric(t(u.initial.first)%*%u.initial.first))
  v.initial.first <- as.numeric(t(k[[i]]) %*% rep(1,p)/p)
  v.initial[[i]] <- v.initial.first /sqrt(as.numeric(t(v.initial.first)%*%v.initial.first))


  # select.paramters.multiple figures out the lambda.u and lambda.v that will maximize correlations
  #  using 
  # We then enter the lambda values into scca.function on the entire dataset.
  # To figure out lambda.u and lambda.v, we run SCCA on training and test data and keep the lambda
  # values which maximize the correlations.  Then we enter the lambda values into scca.function on the
  # entire dataset.  
  parameters <- select.parameters.multiple(X,Y, bound = 1,n.cv = 5, d.vec = d,u.list =  sp.coefs.u.list, v.list =  sp.coefs.v.list, count = i, rob.p = rob)
  lambda.u[i] <- parameters$lambda.u[1]
  lambda.v[i] <- parameters$lambda.v[1]


  uv[[i]] <- scca.function(k[[i]], u.initial[[i]], v.initial[[i]], lambda.u[i], lambda.v[i])

  
  sp.coefs.u.list[[i]]<-uv[[i]]$u.new
  sp.coefs.v.list[[i]]<-uv[[i]]$v.new


  d[i] <- as.numeric(t(sp.coefs.u.list[[i]])%*%k[[i]]%*%sp.coefs.v.list[[i]]) #coefficient for next decomposition

  # Next canonical correlations
  if(rob==T){
    sp.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%sp.coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%sp.coefs.v.list[[i]], method = "spearman")
  }else{
    sp.cor[i] <- cor(X%*%diag(sqrt(1/diag(var.x)))%*%sp.coefs.u.list[[i]],Y%*%diag(sqrt(1/diag(var.y)))%*%sp.coefs.v.list[[i]], method = "pearson")
  }

}  #Close for loop



  
# The actual coefficients given the best lambda values
sp.coef.u <- data.frame(matrix(unlist(sp.coefs.u.list), nrow=length(sp.coefs.u.list[[1]]), byrow=F))
sp.coef.v <- data.frame(matrix(unlist(sp.coefs.v.list), nrow=length(sp.coefs.v.list[[1]]), byrow=F))
  
  

# Turning lists into data frames 

return(list(sp.coef.u = sp.coef.u, sp.coef.v = sp.coef.v, sp.cor = sp.cor, lambda.u = lambda.u, lambda.v = lambda.v))

}
