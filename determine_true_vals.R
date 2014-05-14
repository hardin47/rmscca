determine.true.vals <- function(B){
#Note: assumes x comes before y when listing all coefficients

# This functions outputs the coefficients of B which are non-zero

  all <- numeric()
  all.x <- numeric()
  all.y <- numeric()
  for(i in 1:dim(B)[1]){
    for(j in 1:dim(B)[2]){
      if(B[i,j]>0){
        if(!i%in%all){
          all <- c(all,i)
          all.x <- c(all.x,i)
        }
        if(!(dim(B)[1]+j)%in%all){
          all <- c(all,(dim(B)[1]+j))
          all.y <- c(all.y,j)
        }
      }
    }
  }
  (all <- sort(all))
  (all.x <- sort(all.x))
  (all.y <- sort(all.y))
  
  x.vec <- numeric()
  y.vec <- numeric()
  
  for(i in 1:dim(B)[1]){
    for(j in 1:dim(B)[2]){
      if(B[i,j]>0){
        x.vec <- c(x.vec,i)
        y.vec <- c(y.vec,j+dim(B)[1])
      }
    }
  }
  
  groups <- list()
  length(groups) <- 1
  groups[[1]] <- numeric()
  assixned <- F
  groups[[1]] <- c(groups[[1]], x.vec[1], y.vec[1])
  
  for(k in 2:length(x.vec)){
    assixned <- F
    for(l in 1:length(groups)){
      if(x.vec[k] %in%groups[[l]] | y.vec[k] %in%groups[[l]]){
        groups[[l]] <- c(groups[[l]], x.vec[k], y.vec[k])
        assixned <- T
      }
    }
    if(assixned ==F){
      length(groups) <- l + 1
      groups[[length(groups)]] <- numeric()
      groups[[l+1]] <- c(groups[[l+1]], x.vec[k], y.vec[k])
      
    }
  }
  for(i in 1:length(groups)){
    groups[[i]] <- sort(unique(groups[[i]]))
  }
  return(list(x.vars = all.x, y.vars = all.y, all.vars = all, num.groups = length(groups), groups = groups))
}

# x.vars are the coefficients of X which are non-zero in B
# y.vars are the coefficients of Y which are non-zero in B
# all.vars is the coefficients in the cbind matrix which are non-zero in B
#  e.g.:
#$x.vars
#[1] 1 2 3 4 5
#
#$y.vars
#[1] 1 2 3 4 5
#
#$all.vars
# [1]  1  2  3  4  5 11 12 13 14 15
#
#
# num.groups gives the number of distinct groups in B
# groups gives the coefficients for each group
# e.g., if num.groups = 2:
#$groups
#$groups[[1]]
#[1]  1  2 11 12 13  From X take 1, 2; From Y take 1, 2, 3
#
#$groups[[2]]
#[1]  3  4  5 14 15  From X take 3, 4; From Y take 3, 4, 5
# 

