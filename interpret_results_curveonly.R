interpret.results.curve <- function(output, sp.permcors){

# this function takes as input the output from the results function
# note that we have output & permuted cutoffs from the sparse CCA 

# correlation is the variable which gives the canonical correlations for each pair within the output
# we can subset the results to only consider values which are above a particular cutoff value (output.cor).
# Or we can take the first 5 canonical correlations (which aren't necessarily ranked from 1 to 5), (output.five).

####################
#### SPARSE CCA ####
####################

# In sp.cutoff.val, the correlation for each canonical pair is used.  The permuted corrleations are considered
# to be a null cutoff value for each of the correlations (first, second, third, etc.) within the same dataset.
# Each dataset will have its own cutoff value (as opposed to using one correlation across all datasets).
  sp.cutoff.val = ifelse(sp.permcors[1]<output$sp.Correlation[1], ifelse(sum(sp.permcors < output$sp.Correlation, na.rm=T) == sum(!is.na(output$sp.Correlation)), sum(!is.na(output$sp.Correlation)), min(which(sp.permcors > output$sp.Correlation)) - 1), NA)
  if(!is.na(sp.cutoff.val)) {
	sp.output.perm <- output[1:sp.cutoff.val,]   }else{
	sp.output.perm <- NA}
  
# Num.Not.Trival counts the number of canonical pairs which are not completely zero.
# the which function lists the row numbers, and we take the length of the vector

  sp.Num.Not.Trivial = length(which(output$sp.Num_Total_Pos>0))
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.Num.Not.Trivial.Perm = length(which(sp.output.perm$sp.Num_Total_Pos>0))}else{
  sp.Num.Not.Trivial.Perm = NA}
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
    sp.Num.Can.Pairs.Perm =  sp.cutoff.val # number of CC pairs 
    }else{ sp.Num.Can.Pairs.Perm = 0  }

 
# TP measures the total number of coefficients that are TRUE (with double counting)
# divided by the sum of (total number of coefficients that are TRUE) + (total non-zero coefs that
# aren't in any group)  = proportion of non-zero coefficients which are positive
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.TP.Percent.Perm <- with(sp.output.perm,sum(sp.Num_True_Pos, na.rm=T)/sum(sp.Num_Total_Pos, na.rm=T)) }else{
	sp.TP.Percent.Perm = NA}

	
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.ave.TP.rate.Perm <- with(sp.output.perm,mean(sp.Num_True_Pos/sp.Num_Total_Pos, na.rm=T)) }else{
	sp.ave.TP.rate.Perm = NA}

# Another measure of TP:  sum whether there are (a) ANY true coefficients, and (b) CG
# need to divide by the number of canonical pairs
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.TP.any.Perm <- with(sp.output.perm,sum(sp.Num_True_Pos!=0, na.rm=T)) }else{
	sp.TP.any.Perm = NA}

  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.TP.CG.Perm <- with(sp.output.perm,sum(sp.Complete_Groups!=0, na.rm=T)) }else{
	sp.TP.CG.Perm = NA}


# Seems more consistent with the literature to use FDR = 1-FP...
  sp.FDR.Perm <- 1-sp.TP.Percent.Perm

# Average correlation of the paris we keep  
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.Avg.Cor.Perm <- mean(abs(sp.output.perm$sp.Correlation), na.rm=T)  }else{
	sp.Avg.Cor.Perm = NA  }
  
#  Complete.Groups <- sum(output$Complete_Groups, na.rm=T)
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
    sp.Complete.Groups.Perm <- sum(sp.output.perm$sp.Complete_Groups, na.rm=T) }else{
    sp.Complete.Groups.Perm = NA}
  
# Num_False_Pos measures the number of false positives in a particular canonical pair
# Avg.FP averages the number of false positives over the canonical pairs of interest
# average is taken over the set of canonical pairs
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.Avg.FP.Perm <- mean(sp.output.perm$sp.Num_False_Pos[!is.na(sp.output.perm$sp.Correlation)], na.rm=T)  }else{
	sp.Avg.FP.Perm = NA  }

# Num_True_Pos measures the number of true positives in a particular canonical pair
# Avg.TP averages the number of true positives over the canonical pairs of interest
# average is taken over the set of canonical pairs
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
  sp.Avg.TP.Perm <- mean(sp.output.perm$sp.Num_True_Pos[!is.na(sp.output.perm$sp.Correlation)], na.rm=T)  }else{
	sp.Avg.TP.Perm = NA  }

  
# FN are stored canonical pair by canonical pair.  We track which ones haven't come up yet.
# So, we only need to look for the number of true coefficients that haven't yet been allocated
# at the pair determined by the canonical correlation (5, perm, etc.)
# Need to divide by 48 if B=1, 200 if B=2
  if(!is.na(sp.cutoff.val) & sp.cutoff.val > 1) {
      sp.FN.Perm <- sp.output.perm$sp.Num_False_Neg[dim(sp.output.perm)[1]] }else{
      sp.FN.Perm = NA  }

#################################################################################################

vec <-  numeric(13)
  vec[1] <- sp.Num.Not.Trivial
  vec[2] <- sp.Num.Not.Trivial.Perm
  vec[3] <- sp.TP.Percent.Perm
  vec[4] <- sp.TP.any.Perm
  vec[5] <- sp.TP.CG.Perm
  vec[6] <- sp.ave.TP.rate.Perm
  vec[7] <- sp.Avg.Cor.Perm
  vec[8] <- sp.Complete.Groups.Perm
  vec[9] <- sp.Avg.FP.Perm
  vec[10] <- sp.Avg.TP.Perm
  vec[11] <- sp.FN.Perm
  vec[12] <- sp.FDR.Perm
  vec[13] <- sp.Num.Can.Pairs.Perm  
  return(vec)
}
