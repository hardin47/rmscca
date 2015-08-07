
# to run: R CMD BATCH nullSimScript.R & 

date()

library(foreach)
library(doParallel)
#library(mvnfast)

#setwd("C:/Users/Jo Hardin/Desktop/rsmcca")
source("Final_funcs/build_B.R")
source("Final_funcs/sim_setup.R")
source("Final_funcs/sample_sigma12_function.R")
source("Final_funcs/scca_function.R")
source("Final_funcs/scca_CVperm.R")
source("Final_funcs/Cov_suped.R")

start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")

##  We first set the parameters for running in parallel as well as the 
##  population parameter set-up

num.cluster1 = 25
num.runs1 = 4*num.cluster1

k = 0   # not:, set k=0 to create no relationship between X and Y
p = 500
q = 1000
Btype = 2
num.obs = 50
n.pair = 2  # should be at least 10
nperm=100

# cutoff.perc tells where to cutoff for permutation values
cutoff.perc = 0.9

cor.suped = .5		# the cor of internal X and internal Y
noise = "clean"
# options are clean, t, sym, and asym  (with t, sym, and asym you need noise.level)
# t uses df=2


B <- build.B(k,p,q,Btype)

run1 <- list()
length(run1)<- num.runs1 #need this value

c1 <- makeCluster(num.cluster1)
registerDoParallel(c1)

##  This loop runs the entire simulation in parallel for each dataset.

results.sim <- foreach(i=1:num.runs1, .combine='rbind') %dopar%{

library(mvnfast)

#set.seed(47)
simdata = sim.setup(n.obs = num.obs, B, var.cor=cor.suped, noisetype = noise, Btype=Btype)
sim.output = scca.CVperm(simdata, n.pair, nperm)


# using the permuted correlations to create a curve to determine significance cutoffs
perm.cor.s = sim.output$perm.cor.s
perm.s.curve = apply(perm.cor.s, 2, quantile, probs=cutoff.perc)

perm.cor.p = sim.output$perm.cor.p
perm.p.curve = apply(perm.cor.p, 2, quantile, probs=cutoff.perc)


c(perm.s.curve, sim.output$cor.test.s, perm.p.curve, sim.output$cor.test.p)
}  

fname = paste("nullB",Btype,".n",num.obs,".p",p,".q",q,".",noise,".txt",sep="")
write.table(results.sim, file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

end1 <- date()
end1 <- strptime(end1, "%a %b %d %H:%M:%S %Y")
dif1 <- as.numeric(difftime(end1,start,units="mins"))  # how long the first loop takes, in minutes


write.table(cbind(dif1, dif1, fname), file="times.txt", row.names=F, col.names=F, quote=F, sep="\t", append=T)
