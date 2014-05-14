
results <- function(res, B.mat, bic){
  # res are the results from scca.multiple  

  source("Final_funcs/result_helpers.R")
  source("Final_funcs/determine_true_vals.R")
  
  #Here true is the output of determine.true.vals()
  true <- determine.true.vals(B.mat)
  groups <- true$groups
  full.group = numeric(length(groups))
  #now two options
  if(bic ==T){
    all.coefs <- rbind(abs(rbind(res$bic.coef.u,res$bic.coef.v)))
    cor.vals <- res$bic.cor
    cor.order <- order(abs(res$bic.cor),decreasing =T)
  }else{
    all.coefs <- rbind(abs(rbind(res$nobic.coef.u,res$nobic.coef.v)))
    cor.vals <- res$nobic.cor
    cor.order <- order(abs(res$nobic.cor),decreasing =T)
  }
  
  num.pairs = dim(all.coefs)[2]
  num.vars = dim(all.coefs)[1]
  
  
  nonzero = list()
  length(nonzero) = num.pairs
  for(i in 1:num.pairs){
    nonzero[[i]] = which(all.coefs[,i]>.00001)
  }

  
  
  
  #Univariate TP - # true vars represented in all pairs, with double counting
  #Univariate FP - # non-zero coefs that aren't in any group
  #Univariate FN - # true vals that never came up
  
  #Creates a summary table for each canonical pair
  big <- data.frame(nrow = num.pairs, ncol = 9)
  true.vals <- c()
  for(i in 1:num.pairs){
    UTP = length(which(true$all%in%nonzero[[cor.order[i]]]))
    UFP = length(which(!nonzero[[cor.order[i]]]%in%true$all))
    pos <- true$all[which(true$all.vars%in%nonzero[[cor.order[i]]])]
    true.vals <- c(true.vals, pos)
    true.vals <- sort(unique(true.vals))
    UFN <- length(true$all) - length(true.vals)
    big[i,1] = complete.groups(groups,nonzero[[cor.order[i]]])
    big[i,2] = UTP
    big[i,3] = UFP
    big[i,4] = UTP+UFP
    big[i,5] = UTP/UFP
    big[i,6] = paste(pos,collapse = ",")
    big[i,7] = UFN
    big[i,8] = round(cor.vals[cor.order[i]],3)
    big[i,9] = cor.order[i]
  }
  for(i in 1:num.pairs){
    
  }

  names(big) <- c("Complete_Groups", "Num_True_Pos", "Num_False_Pos","Num_Total_Pos", "Ratio_TP/FP", "True_Pos_Values","Num_False_Neg", "Correlation", "Cor.Order")
#  big <- big[order(abs(big$Correlation), decreasing = T),]

    
  return(big)
}
