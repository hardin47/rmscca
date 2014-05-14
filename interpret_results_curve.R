interpret.results.curve <- function(output, cor.val, permcors){

# this function takes as input the output from the results function

  output.cor <- subset(output, Correlation>cor.val)
  output.five <- output[1:5,]

  cutoff.val = ifelse(permcors[1]<output$Correlation[1], ifelse(sum(permcors < output$Correlation, na.rm=T) == sum(!is.na(output$Correlation)), sum(!is.na(output$Correlation)), min(which(permcors > output$Correlation)) - 1), NA)
  if(!is.na(cutoff.val)) {
	output.perm <- output[1:cutoff.val,]   }else{
	output.perm <- NA}
  
# Num.Not.Trival counts the number canonical pairs which are not completely zero.
# the which lists the row numbers, and we take the length of the vector

  Num.Not.Trivial = length(which(output$Num_Total_Pos>0))
  Num.Not.Trivial.Perm = 
  Num.Can.Pairs.Perm = cutoff.val

 
# TP measures the total number of coefficients that are TRUE (with double counting)
# divided by the sum of (total number of coefficients that are TRUE) + (total non-zero coefs that
# aren't in any group)  = proportion of non-zero coefficients which are positive
  TP.Percent.Tot <- with(output,sum(Num_True_Pos, na.rm=T)/sum(Num_Total_Pos, na.rm=T))
  TP.Percent.Five <- with(output.five,sum(Num_True_Pos, na.rm=T)/sum(Num_Total_Pos, na.rm=T))
  TP.Percent.Cor <- with(output.cor,sum(Num_True_Pos, na.rm=T)/sum(Num_Total_Pos, na.rm=T))
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  TP.Percent.Perm <- with(output.perm,sum(Num_True_Pos, na.rm=T)/sum(Num_Total_Pos, na.rm=T)) }else{
	TP.Percent.Perm = NA}


# Seems more consistent with the literature to use FDR = 1-FP...
  FDR.Tot <- 1-TP.Percent.Tot  
  FDR.Five <- 1-TP.Percent.Five
  FDR.Cor <- 1-TP.Percent.Cor 
  FDR.Perm <- 1-TP.Percent.Perm

  
  Avg.Cor.Tot <- mean(abs(output$Correlation), na.rm=T)
  Avg.Cor.Five <- mean(abs(output.five$Correlation), na.rm=T)
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  Avg.Cor.Perm <- mean(abs(output.perm$Correlation), na.rm=T)  }else{
	Avg.Cor.Perm = NA  }
  Num.vec.Cor <- dim(output.cor)[1]
  
  Complete.Groups <- sum(output$Complete_Groups, na.rm=T)
  
# Num_False_Pos measures the number of false positives in a particular canonical pair
# Avg.FP averages the number of false positives over the canonical pairs of interest
  Avg.FP.Tot <- mean(output$Num_False_Pos[!is.na(output$Correlation)], na.rm=T)
  Avg.FP.Five <- mean(output.five$Num_False_Pos[!is.na(output.five$Correlation)], na.rm=T)
  Avg.FP.Cor <- mean(output.cor$Num_False_Pos[!is.na(output.cor$Correlation)], na.rm=T)
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  Avg.FP.Perm <- mean(output.perm$Num_False_Pos[!is.na(output.perm$Correlation)], na.rm=T)  }else{
	Avg.FP.Perm = NA  }


  
# FN are stored canonical pair by canonical pair.  We track which ones haven't come up yet.
# So, we only need to look for the number of true coefficients that haven't yet been allocated
# at the pair determined by the canonical correlation (5, perm, etc.)
  FN.Tot <- output$Num_False_Neg[dim(output)[1]]
  FN.Five <- output$Num_False_Neg[5]
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  FN.Perm <- output$Num_False_Neg[cutoff.val]  }else{
	FN.Perm = NA  }
  if(dim(output.cor)[1]>0){
    FN.Cor <-output.cor$Num_False_Neg[dim(output.cor)[1]]
  }else{
    FN.Cor <- NA
  }
  vec <-  numeric(23)
  vec[1] <- Num.Not.Trivial
  Vec[2] <- Num.Not.Trivial.Perm
%  vec[2] <- TP.Percent.Tot
%  vec[3] <- TP.Percent.Five
%  vec[4] <- TP.Percent.Cor
  vec[3] <- TP.Percent.Perm
%  vec[6] <- Avg.Cor.Tot
%  vec[7] <- Avg.Cor.Five
%  vec[8] <- Num.vec.Cor
  vec[4] <- Avg.Cor.Perm
%  vec[10] <- Complete.Groups
  vec[5] <- Compelte.Groups.Perm
%  vec[11] <-  Avg.FP.Tot
%  vec[12] <- Avg.FP.Five
%  vec[13] <- Avg.FP.Cor
  vec[6] <- Avg.FP.Perm
%  vec[15] <-  FN.Tot
%  vec[16] <-  FN.Five
%  vec[17] <-  FN.Cor
  vec[7] <-  FN.Perm
%  vec[19] <- FDR.Tot
%  vec[20] <- FDR.Five
%  vec[21] <- FDR.Cor
  vec[8] <- FDR.Perm
  vec[9] <- Num.Can.Pairs.Perm  
  return(vec)
}
