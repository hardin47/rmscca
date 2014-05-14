interpret.results <- function(output, cor.val){
  output.cor <- subset(output, Correlation>cor.val)
  output.five <- output[1:5,]
  
  Num.Not.Trivial = max(which(output$Num_Total_Pos>0))
  
  TP.Percent.Tot <- with(output,sum(Num_True_Pos)/sum(Num_Total_Pos))
  TP.Percent.Five <- with(output.five,sum(Num_True_Pos)/sum(Num_Total_Pos))
  TP.Percent.Cor <- with(output.cor,sum(Num_True_Pos)/sum(Num_Total_Pos))
  
  Avg.Cor.Tot <- mean(abs(output$Correlation))
  Avg.Cor.Five <- mean(abs(output.five$Correlation))
  Num.vec.Cor <- dim(output.cor)[1]
  
  Complete.Groups <- sum(output$Complete_Groups)
  
  Avg.FP.Tot <- mean(output$Num_False_Pos)
  Avg.FP.Five <- mean(output.five$Num_False_Pos)
  Avg.FP.Cor <- mean(output.cor$Num_False_Pos)
  
  FN.Tot <- output$Num_False_Neg[dim(output)[1]]
  FN.Five <- output$Num_False_Neg[5]
  if(dim(output.cor)[1]>0){
    FN.Cor <-output.cor$Num_False_Neg[dim(output.cor)[1]]
  }
  else{
    FN.Cor <- NA
  }
  vec <-  numeric(12)
  vec[1] <- Num.Not.Trivial
  vec[2] <- TP.Percent.Tot
  vec[3] <- TP.Percent.Five
  vec[4] <- TP.Percent.Cor
  vec[5] <- Avg.Cor.Tot
  vec[6] <- Avg.Cor.Five
  vec[7] <- Num.vec.Cor
  vec[8] <- Complete.Groups
  vec[9] <-  Avg.FP.Tot
  vec[10] <- Avg.FP.Five
  vec[11] <- Avg.FP.Cor
  vec[12] <-  FN.Tot
  vec[13] <-  FN.Five
  vec[14] <-  FN.Cor
  
  return(vec)
}
