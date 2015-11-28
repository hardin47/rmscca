complete.groups <- function(groups,pos){
  # Num of complete groups
  full.group <-0
  for(j in 1:length(groups)){
    if(sum(groups[[j]]%in%pos)==length(groups[[j]])){
      full.group = full.group + 1    
    }
  }
  return(full.group)
  
}

# this function is used in results.R to give the first column of the big matrix
#
# complete.groups(groups, nonzero[[cor.order[i]]]) for the ith canonical pair
#
# cor.order[i] is the order of the canonical pair based on the sorted can corr
#
#    nonzero[[i]] = which(all.coefs[,i]>.00001)
#
# so - for each canonical vector - the complete.groups function takes the intersection of 
# the non-zero coefficients with the predefined groups.  If an entire group is contained
# in the non-zero coefficients, then complete groups gets incremented by one


