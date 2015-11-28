## build.B returns a pxq diagonal block matrix. Block elements are k or 0. ##
## B defines the relationship between X and Y. E(Y) = B*X		   ##


build.B <- function(k,p,q, Btype){
# Called by fullSimScript.R and nullSimScript.R
# Calls no other rmscca functions

B <- matrix(0,nrow = p, ncol = q)	#Create the initial pxq matrix, all elements 0.


# Btype = 0: 2x3, 3x2

if(Btype==0){
  for(i in 1:2){
    for(j in 1:3){
      B[i,j] = k	#Rows 1-2, Collumns 1-3, entries of B assigned k.
    }
  }
  for(i in 3:5){
    for(j in 4:5){
      B[i,j] = k	#Rows 3-5, Columns 4-5, entries of B assigned k.
    }
  }
}

# Btype = 1: 12x8, 8x4, 12x4

if (Btype==1){
  for(i in 1:12){
    for(j in 1:8){
      B[i,j] <- k	#Rows 1-12, Columns 1-8, entries of B assigned k.
    }
  }
  for(i in 13:20){
    for(j in 9:12){
      B[i,j] <- k	#Rows 13-20, Columns 9-12, " ".
    }
  }
  for(i in 21:32){
    for(j in 13:16){
      B[i,j] <- k	#Rows 21-32, Columns 13-16, " ".
    }
  }
}

# Btype = 2: 10x20, 5x5, 20x10, 50x50, 15x15

if (Btype==2){
  for(i in 1:10){
    for(j in 1:20){
      B[i,j] <- k	#Rows 1-10, Columns 1-20, entries of B assigned k.
    }
  }
  for(i in 11:16){
    for(j in 21:26){
      B[i,j] <- k	#Rows 11-16, Columns 21:26, " ".
    }
  }
  for(i in 17:36){
    for(j in 27:36){
      B[i,j] <- k	#Rows 17-36, Columns 27-36, " ".
    }
  }
  for(i in 37:86){
    for(j in 37:86){
      B[i,j] <- k	#Rows 37-86, Columns 37-86, " ".
    }
  }
  for(i in 87:100){
    for(j in 87:100){
      B[i,j] <- k	#Rows 87-100, Columns 87-100, " ".
    }
  }
}


return(B)

}