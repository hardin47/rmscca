date()

library(foreach)
library(doParallel)

source("Final_funcs/results.R")
source("Final_funcs/interpret_results.R")
source("Final_funcs/interpret_results_curve.R")
source("Final_funcs/build_B.R")
source("Final_funcs/determine_true_vals.R")
source("Final_funcs/sim_setup_noise_suped.R")
#source("Final_funcs/sim_setup_noise.R")
#source("Final_funcs/sim_setup_suped.R")
#source("Final_funcs/sim_setup.R")
source("Final_funcs/Hardin_Rand_Cov_Func.R")
source("Final_funcs/randCov_suped.R")

source("Final_funcs/sample.sigma12.function.R")
#source("Final_funcs/sample.sigma12.function.fullsig.R")
source("Final_funcs/select.parameters.multiple.R")
source("Final_funcs/bic.function.R")
source("Final_funcs/scca.function.R")
source("Final_funcs/scca.multiple.R")


start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")

##  We first set the parameters for running in parallel as well as the 
##   population parameter set-up

num.cluster1 = 25
num.runs1 = 4*num.cluster1

num.cluster2 = 25
num.runs2 = 4*num.cluster2

k = 1
p = 100
q = 100
Btype = 1
num.obs = 100

cor.suped = .3
noise.level = 0.1
noise = "sym"
# Options are clean, t, sym, and asym  (with sym and asym you need noise.level)
# t uses df=3, we might want a lower df? 1?


B = build.B(k,p,q,Btype)

# cutoff parameters for later
cutoff.perc = 0.99
cor.val.cutoff = 0.4

setup = sim.setup.noise.suped(num.obs, B, contamination=noise.level, var.cor=cor.suped, noisetype = noise, Btype=Btype)
X = setup$X
Y = setup$Y

run1 <- list()
length(run1)<- num.runs1 #need this val

c1 <- makeCluster(num.cluster1)
registerDoParallel(c1)

##  This first foreach loop runs one dataset many times (permuted) to give us
##  an estimate for the null cutoff curve.  We use that cutoff curve in the
##  next foreach loop so that we can make a determination of whether the 
##  cannonical pair should be considered as significant.  

perm.results.sim <- foreach(i=1:num.runs1, .combine='cbind') %dopar%{

  Yperm = Y[sample(1:num.obs, replace=F),]

  res.s <- scca.multiple(X, Yperm, rob = T, bic=F)#, fullsig=F) #Here bic=T means compute bic as well as nobic
  res.p <- scca.multiple(X, Yperm, rob = F, bic=F)#, fullsig=F)

#  c(res.s$nobic.cor, res.s$bic.cor, res.p$nobic.cor, res.p$bic.cor)
   c(res.s$nobic.cor, res.p$nobic.cor)

}  

fname = paste("permB",Btype,".n",num.obs,".p",p,".q",q,".",noise,".txt",sep="")
write.table(perm.results.sim, file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

end1 <- date()
end1 <- strptime(end1, "%a %b %d %H:%M:%S %Y")
dif1 <- as.numeric(difftime(end1,start,units="mins"))  # how long the first loop takes, in minutes


##  Before running the second loop, we need to parse out the information from the permutation
##  test to create the cutoff curves


mpq = min(p,q)
perm.s = perm.results.sim[1:mpq,]
perm.s = apply(perm.s, 2, sort, decreasing=T, na.last=T)
perm.s.curve = apply(perm.s,1, quantile, probs=cutoff.perc, na.rm=T)

perm.p = perm.results.sim[(mpq+1):(2*mpq),]
perm.p = apply(perm.p, 2, sort, decreasing=T, na.last=T)
perm.p.curve = apply(perm.p,1, quantile, probs=cutoff.perc, na.rm=T)


##  Now, given the cutoff curves, we can run the full simulation.  Note that we need to generate
##  a new dataset in each loop, but we are always comparing it to the curves above.
##  For each dataset, we compute the number of false positives, false negatives, total group
##  relationships, etc.

###  FIX ABOVE TO REMIND OURSELVES EXACTLY WHAT WE ARE OUTPUTTING


run2 <- list()
length(run2)<- num.runs2 #need this val

c2 <- makeCluster(num.cluster2)
registerDoParallel(c2)

sim.results <- foreach(i=1:num.runs2, .combine='rbind') %dopar%{

  setup <- sim.setup.noise.suped(num.obs,B,contamination = noise.level, var.cor  = cor.suped, noisetype=noise, Btype=Btype)
  X <- setup$X
  Y <- setup$Y
  res.s <- scca.multiple(X,Y, rob = T, bic=F) #Here bic=T means compute bic as well as nobic
  res.p <- scca.multiple(X,Y, rob = F, bic=F)
  

  output.s.nobic <- results(res.s, B, bic=F)
  output.p.nobic <- results(res.p, B, bic = F)


# we are no longer using BIC which should save computational time  
#  output.s.bic <- results(res.s, B, bic=T)
#  output.p.bic <- results(res.p, B, bic = T)
  
#   big.s.bic <- big.s.nobic <- big.p.bic <- big.p.nobic <- matrix(nrow = num.runs2, ncol = 12)
#   
#   big.s.bic[i,] <-  interpret.results(output.s.bic,cor.val.cutoff)
#   big.s.nobic[i,] <-interpret.results(output.s.nobic,cor.val.cutoff)
#   big.p.bic[i,] <- interpret.results(output.p.bic,cor.val.cutoff)
#   big.p.nobic[i,] <- interpret.results(output.p.nobic, cor.val.cutoff)
  
 c( interpret.results.curve(output.s.nobic,cor.val.cutoff,perm.s.curve ),
       interpret.results.curve(output.p.nobic, cor.val.cutoff,perm.p.curve))
}  
    

fname2 = paste("simB",Btype,".n",num.obs,".p",p,".q",q,".",noise,".txt",sep="")
write.table(sim.results, file=fname2, row.names=F, quote=F, col.names=F, sep="\t")

  end2 <- date()
  end2 <- strptime(end2, "%a %b %d %H:%M:%S %Y")
  dif2 <- as.numeric(difftime(end2, end1, units = "mins"))

write.table(c(dif1, dif2, fname2), file="times.txt", row.names=F, col.names=F, quote=F, sep="\t", append=T)