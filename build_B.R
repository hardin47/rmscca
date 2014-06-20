## build.B returns a pxq diagonal block matrix. Block elements are k or 0. ##
## B defines the relationship between X and Y. E(Y) = B*X		   ##


build.B <- function(k,p,q, Btype){
# Called by big_sim_cutoff.R
# Calls no other rmscca functions

B <- matrix(0,nrow = p, ncol = q)	#Create the initial pxq matrix, all elements 0.

if(Btype==0){
  for(i in 1:2){
    for(j in 1:3){
      B[i,j] = k	#Rows 1-2, Collumns 1-3, entries of B assigned k.
    }
  }
  for(i in 3:5){
    for(j in 4:5){
      B[i,j] = k	#Rows 3-5, Collumns 4-5, entries of B assigned k.
    }
  }
}

if (Btype==1){
  for(i in 1:12){
    for(j in 1:8){
      B[i,j] <- k	#Rows 1-12, Collumns 1-8, entries of B assigned k.
    }
  }
  for(i in 13:20){
    for(j in 9:12){
      B[i,j] <- k	#Rows 13-20, Collumns 9-12, " ".
    }
  }
  for(i in 21:32){
    for(j in 13:16){
      B[i,j] <- k	#Rows 21-32, Collumns 13-16, " ".
    }
  }
}

# For Btype == 2 we need to think about a bigger matrix to test.

if(Btype==2)  {}


return(B)

}