scca.function <- function(k, u.initial, v.initial, lambda.u, lambda.v,gamma = 1)
{
  i <- 0    # number of iterations used by SCCA
  eps <- 0.001	# convergence criterion 
  max.iter <- 50	# maximum nuber of iterations
  diff.u <- eps*10	
  diff.v <- eps*10
 # svdvec <- svd(k,nu = 1, nv = 1)
  
  while ((i < max.iter) & ((diff.u > eps) || (diff.v > eps)) )
  {
    i <- i+1
    
    
    
    # Update left singular vector
    
    vx <-  k %*% v.initial
    length.vx <- as.numeric(sqrt(t(vx)%*%vx))
    if(length.vx==0) length.vx <- 1
    vx <- vx / length.vx
 #   svdu <- norm(abs(as.matrix(svdvec$u[,1])^gamma),"F")
    u.new <- abs(vx) - 0.5*lambda.u#/svdu
    u.new <- (u.new + abs(u.new))/2
    u.new <- u.new*sign(vx)
    length.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
    if(length.u.new==0) length.u.new <- 1
    u.new <- u.new / length.u.new
    
    
    # Update right singular vector
    
    ux <- t(k) %*% u.new
    length.ux <- as.numeric(sqrt(t(ux)%*%ux))
    if(length.ux==0) length.ux <- 1
    ux <- ux / length.ux
 #   svdv <- norm(abs(as.matrix(svdvec$v[,1])^gamma),"F")
    v.new <- abs(ux) - 0.5*lambda.v#/svdv
    v.new <- (v.new + abs(v.new))/2
    v.new <- v.new * sign(ux)
    length.v.new <- as.numeric(sqrt(t(v.new)%*%v.new))
    if(length.v.new==0) length.v.new <- 1
    v.new <- v.new / length.v.new
    
    
    # Convergence measures
    
    diff.v <- max(abs(v.initial - v.new))
    diff.u <- max(abs(u.initial - u.new))
    
    v.initial <- v.new
    u.initial <- u.new
  }
  
  # Report the results:
  # u.new is computed left singular vector
  # v.new is computed right singular vector
  # i is the number of iterations used by SCCA
  
  list(u.new=u.new, v.new=v.new, i=i)
}
