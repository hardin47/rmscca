## Runs SCCA for one canonical pair ##

scca.function <- function(k, u.initial, v.initial, lambda.u, lambda.v) {

# Called by scca.multiple.R & select.multiple.parameters.R & scca_CVperm
# Calls no other rmscca functions

i <- 0 		  # Number of iterations used by SCCA
eps <- 0.01	  # Convergence criterion (was 0.001)
max.iter <- 20	  # Maximum nuber of iterations (was 50)
diff.u <- eps*10	
diff.v <- eps*10

while ((i < max.iter) & ((diff.u > eps) || (diff.v > eps)) ) {

  i <- i+1
   
  # Update left singular vector
  vx <-  k %*% v.initial

  length.vx <- as.numeric(sqrt(t(vx)%*%vx))		# Normalize vx
  if(length.vx==0) length.vx <- 1			# 
  vx <- vx / length.vx					#


  u.new <- abs(vx) - 0.5*lambda.u        		# Apply soft-thresholding to
  u.new <- (u.new + abs(u.new))/2			# obtain u.new, a
  u.new <- u.new*sign(vx)				# sparse solution


  length.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))	# Normalize u.new
  if(length.u.new==0) length.u.new <- 1			#
  u.new <- u.new / length.u.new				#
    
    
  # Update right singular vector
  ux <- t(k) %*% u.new

  length.ux <- as.numeric(sqrt(t(ux)%*%ux))		# Normalize ux
  if(length.ux==0) length.ux <- 1			#  
  ux <- ux / length.ux					#


  v.new <- abs(ux) - 0.5*lambda.v       		# Apply soft-thresholding to
  v.new <- (v.new + abs(v.new))/2			# obtain v.new, a
  v.new <- v.new * sign(ux)				# sparse solution


  length.v.new <- as.numeric(sqrt(t(v.new)%*%v.new))	# Normalize v.new
  if(length.v.new==0) length.v.new <- 1			#
  v.new <- v.new / length.v.new				#
    
    
  # Convergence measures  
  diff.v <- max(abs(v.initial - v.new))			# Update diff.v
  diff.u <- max(abs(u.initial - u.new))			# Update diff.u
    
  v.initial <- v.new					# Update v.initial
  u.initial <- u.new					# Update u.initial

  } # Close while loop
  
  # Report the results:
  # u.new is the (sparse) canonical vector corresponding to X variables
  # v.new is the (sparse) canonical vector corresponding to Y variables
  # i is the number of iterations used by SCCA
  
  list(u.new=u.new, v.new=v.new, i=i)
}
