interpret.results.curve <- function(output, cor.val, permcors){

# this function takes as input the output from the results function

# correlation is the variable which gives the canonical correlations for each pair within the output
# we can subset the results to only consider values which are above a particular cutoff value (output.cor).
# Or we can take the first 5 canonical correlations (which aren't necessarily ranked from 1 to 5), (output.five).


# In cutoff.val, the correlation for each canonical pair is used.  The permuted corrleations are considered
# to be a null cutoff value for each of the correlations (first, second, third, etc.) within the same dataset.
# Each dataset will have its own cutoff value (as opposed to using one correlation across all datasets).
  cutoff.val = ifelse(permcors[1]<output$Correlation[1], ifelse(sum(permcors < output$Correlation, na.rm=T) == sum(!is.na(output$Correlation)), sum(!is.na(output$Correlation)), min(which(permcors > output$Correlation)) - 1), NA)
  if(!is.na(cutoff.val)) {
	output.perm <- output[1:cutoff.val,]   }else{
	output.perm <- NA}
  
# Num.Not.Trival counts the number of canonical pairs which are not completely zero.
# the which function lists the row numbers, and we take the length of the vector

  Num.Not.Trivial = length(which(output$Num_Total_Pos>0))
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  Num.Not.Trivial.Perm = length(which(output.perm$Num_Total_Pos>0))}else{
  Num.Not.Trivial.Perm = NA}
  if(!is.na(cutoff.val) & cutoff.val > 1) {
    Num.Can.Pairs.Perm =  cutoff.val # number of CC pairs 
    }else{ Num.Can.Pairs.Perm = 0  }

 
# TP measures the total number of coefficients that are TRUE (with double counting)
# divided by the sum of (total number of coefficients that are TRUE) + (total non-zero coefs that
# aren't in any group)  = proportion of non-zero coefficients which are positive
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  TP.Percent.Perm <- with(output.perm,sum(Num_True_Pos, na.rm=T)/sum(Num_Total_Pos, na.rm=T)) }else{
	TP.Percent.Perm = NA}


# Seems more consistent with the literature to use FDR = 1-FP...
  FDR.Perm <- 1-TP.Percent.Perm

# Average correlation of the paris we keep  
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  Avg.Cor.Perm <- mean(abs(output.perm$Correlation), na.rm=T)  }else{
	Avg.Cor.Perm = NA  }
  
#  Complete.Groups <- sum(output$Complete_Groups, na.rm=T)
  if(!is.na(cutoff.val) & cutoff.val > 1) {
    Complete.Groups.Perm <- sum(output.perm$Complete_Groups, na.rm=T) }else{
    Complete.Groups.Perm = NA}
  
# Num_False_Pos measures the number of false positives in a particular canonical pair
# Avg.FP averages the number of false positives over the canonical pairs of interest
  if(!is.na(cutoff.val) & cutoff.val > 1) {
  Avg.FP.Perm <- mean(output.perm$Num_False_Pos[!is.na(output.perm$Correlation)], na.rm=T)  }else{
	Avg.FP.Perm = NA  }


  
# FN are stored canonical pair by canonical pair.  We track which ones haven't come up yet.
# So, we only need to look for the number of true coefficients that haven't yet been allocated
# at the pair determined by the canonical correlation (5, perm, etc.)
  if(!is.na(cutoff.val) & cutoff.val > 1) {
      FN.Perm <- output.perm$Num_False_Neg[dim(output.perm)[1]] }else{
      FN.Perm = NA  }

vec <-  numeric(9)
  vec[1] <- Num.Not.Trivial
  vec[2] <- Num.Not.Trivial.Perm
  vec[3] <- TP.Percent.Perm
  vec[4] <- Avg.Cor.Perm
  vec[5] <- Complete.Groups.Perm
  vec[6] <- Avg.FP.Perm
  vec[7] <-  FN.Perm
  vec[8] <- FDR.Perm
  vec[9] <- Num.Can.Pairs.Perm  
  return(vec)
}
