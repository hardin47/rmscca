build.B <- function(k,p,q, Btype){
  B <- matrix(rep(0,p*q),nrow = p, ncol = q)
  
  if(Btype==0){
    for(i in 1:2){
      for(j in 1:3){
        B[i,j] = k
      }
    }
    for(i in 3:5){
      for(j in 4:5){
        B[i,j] = k
      }
    }
    
    
  } 
    
    if (Btype==1){
      for(i in 1:12){
        for(j in 1:8){
          B[i,j] <- k
        }
      }
      for(i in 13:20){
        for(j in 9:12){
          B[i,j] <- k
        }
      }
      for(i in 21:32){
        for(j in 13:16){
          B[i,j] <- k
        }
      }
    }
    
    if(Btype==2)  {}
    
    return(B)
  }