
results <- function(res, B.mat, n.pair){
  # res are the results from scca.multiple  

  source("Final_funcs/result_helpers.R")
  source("Final_funcs/determine_true_vals.R")
  
  #Here true is the output of determine.true.vals()
  true <- determine.true.vals(B.mat)
  groups <- true$groups
  full.group = numeric(length(groups))
#  num.pairs = min(dim(B.mat))
  num.pairs = n.pair

  # get results from SCCA and CCA
    sp.coefs <- rbind(abs(rbind(res$sp.coef.u,res$sp.coef.v)))
    sp.cor.vals <- res$sp.cor
    sp.cor.order <- order(abs(res$sp.cor),decreasing =T)
  
  sp.nonzero = list()
  length(sp.nonzero) = num.pairs
  for(i in 1:num.pairs){
    sp.nonzero[[i]] = which(abs(sp.coefs[,i])>.00001)
  }
  
  
  
  #Univariate TP - # true vars represented in all pairs, with double counting
  #Univariate FP - # non-zero coefs that aren't in any group
  #Univariate FN - # true vals that never came up
  
  #Creates a summary table for each canonical pair
  big <- data.frame(nrow = num.pairs, ncol = 9)
  sp.true.vals <- c()
  for(i in 1:num.pairs){
    sp.UTP = length(which(true$all%in%sp.nonzero[[sp.cor.order[i]]]))
    sp.UFP = length(which(!sp.nonzero[[sp.cor.order[i]]]%in%true$all))
    sp.pos <- true$all[which(true$all.vars%in%sp.nonzero[[sp.cor.order[i]]])]
    sp.true.vals <- sort(unique(c(sp.true.vals, sp.pos)))
    sp.UFN <- length(true$all) - length(sp.true.vals)
    big[i,1] = complete.groups(groups,sp.nonzero[[sp.cor.order[i]]])
    big[i,2] = sp.UTP
    big[i,3] = sp.UFP
    big[i,4] = sp.UTP+sp.UFP
    big[i,5] = sp.UTP/sp.UFP
    big[i,6] = paste(sp.pos,collapse = ",")
    big[i,7] = sp.UFN
    big[i,8] = round(sp.cor.vals[sp.cor.order[i]],3)
    big[i,9] = sp.cor.order[i]
  }


  names(big) <- c("sp.Complete_Groups","sp.Num_True_Pos", "sp.Num_False_Pos","sp.Num_Total_Pos", "sp.Ratio_TP/FP", "sp.True_Pos_Values","sp.Num_False_Neg", "sp.Correlation", "sp.Cor.Order")
#  big <- big[order(abs(big$Correlation), decreasing = T),]

    
  return(big)
}
