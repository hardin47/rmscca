
# to run: R CMD BATCH fullSimScript.R & 

date()

library(foreach)
library(doParallel)


source("Final_funcs/build_B.R")
source("Final_funcs/sim_setup.R")
source("Final_funcs/Cov_suped.R")
source("Final_funcs/sample_sigma12_function.R")
source("Final_funcs/scca_function.R")
source("Final_funcs/scca_CVperm.R")

source("Final_funcs/result_helpers.R")
source("Final_funcs/determine_true_vals.R")
source("Final_funcs/results.R")
source("Final_funcs/interpret_results_curveonly.R")

start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")

##  We first set the parameters for running in parallel as well as the 
##  population parameter set-up

num.cluster1 = 25
num.runs1 = 4*num.cluster1

k = 1
p = 500
q = 1000
Btype = 2
num.obs = 500
n.pair = 10  # should be at least 10
nperm=100

# cutoff.perc tells where to cutoff for permutation values
cutoff.perc = 0.9

cor.suped = .2		# the cor of internal X and internal Y
noise = "clean"
# options are clean, t, sym, and asym  (with t, sym, and asym you need noise.level)
# t uses df=2, we might want a lower df? 1?


B <- build.B(k,p,q,Btype)
#set.seed(47)

run1 <- list()
length(run1)<- num.runs1 #need this val

c1 <- makeCluster(num.cluster1)
registerDoParallel(c1)

##  This loop runs the entire simulation in parallel for each dataset.

results.sim <- foreach(i=1:num.runs1, .combine='rbind') %dopar%{

library(mvnfast)

simdata = sim.setup(num.obs, B, var.cor=cor.suped, noisetype = noise, Btype=Btype)
sim.output = scca.CVperm(simdata, n.pair, nperm)


# using the permuted correlations to create a curve to determine significance cutoffs
perm.cor.s = sim.output$perm.cor.s
perm.s.curve = apply(perm.cor.s, 2, quantile, probs=cutoff.perc)

perm.cor.p = sim.output$perm.cor.p
perm.p.curve = apply(perm.cor.p, 2, quantile, probs=cutoff.perc)


# mapping new output to the same form as the previous output
res.s = list()
res.s$sp.coef.u = data.frame(matrix(unlist(sim.output$alphas.s), nrow=length(sim.output$alphas.s[[1]]), byrow=F))
res.s$sp.coef.v = data.frame(matrix(unlist(sim.output$betas.s), nrow=length(sim.output$betas.s[[1]]), byrow=F))
res.s$sp.cor = sim.output$cor.test.s

res.p = list()
res.p$sp.coef.u = data.frame(matrix(unlist(sim.output$alphas.p), nrow=length(sim.output$alphas.p[[1]]), byrow=F))
res.p$sp.coef.v = data.frame(matrix(unlist(sim.output$betas.p), nrow=length(sim.output$betas.p[[1]]), byrow=F))
res.p$sp.cor = sim.output$cor.test.p


# counting false positives, false negatives, etc.
output.s <- results(res.s, B, n.pair)
output.p <- results(res.p, B, n.pair)

c( interpret.results.curve(output.s, perm.s.curve ),
       interpret.results.curve(output.p, perm.p.curve), sim.output$lambda1.s, sim.output$lambda1.p)

}  

fname = paste("secondsimB",Btype,".n",num.obs,".p",p,".q",q,".",noise,".txt",sep="")
write.table(results.sim, file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

end1 <- date()
end1 <- strptime(end1, "%a %b %d %H:%M:%S %Y")
dif1 <- as.numeric(difftime(end1,start,units="mins"))  # how long the first loop takes, in minutes


write.table(cbind(dif1, dif1, fname), file="times.txt", row.names=F, col.names=F, quote=F, sep="\t", append=T)
