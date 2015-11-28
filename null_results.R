
n.obs=c(50,100,500,1000)
pq = rbind(c(100,100),c(500,1000))
noise=c("clean","sym","t","asym", "noise")
which.B = c("B1","B2")
fullnull.sim=c()

for(i in c(1:4)){  # num observations
  for (j in 2){ # pq combo
    for(k in c(1)){  # noise
      for(l in 2){  # B type
        
        filename = paste("~/Desktop/rsmcca/Final_Sims/null",which.B[l],".n",n.obs[i],".p",pq[j,1],".q",pq[j,2],".",noise[k],".txt", sep="")
        
        # null datat is 100x8;  100 reps; first 4 RMSCCA & second 4 MSCCA
        # first pair is 2 CCs permutation curve, next 2 are actual correlations
        # type I error is if [1,3] > [1,1]
        null.new = read.delim(filename, header=FALSE)
        names(null.new) = c("scurve1", "scurve2", "scor1", "scor2",
                            "pcurve1", "pcurve2", "pcor1", "pcor2")
        reps = dim(null.new)[1]
        
        # full results
        
        stypeI <- sum(null.new$scurve1 < null.new$scor1) / reps
        ptypeI <- sum(null.new$pcurve1 < null.new$pcor1) / reps
        
        null.typeI = cbind(stypeI, ptypeI, "n.obs" = n.obs[i], "p" = pq[j,1], "q" = pq[j,2], "noise"=noise[k], "Btype"=which.B[l])
        
        fullnull.sim = rbind(fullnull.sim, null.typeI)
        
      } }}}

fullnull.sim


