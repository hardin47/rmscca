
# to run: R CMD BATCH running_scca_CVperm.R & 

date()

library(foreach)
library(doParallel)
library(Matrix)
library(MASS)
library(mixstock)


source("Final_funcs/results.R")
source("Final_funcs/interpret_results_curveonly.R")
source("Final_funcs/build_B.R")
source("Final_funcs/determine_true_vals.R")
source("Final_funcs/sim_setup_noise_suped.R")
source("Final_funcs/result_helpers.R")
source("Final_funcs/determine_true_vals.R")
source("Final_funcs/sample.sigma12.function.R")
source("Final_funcs/select.parameters.multiple.R")
source("Final_funcs/scca.function.R")
source("Final_funcs/scca.multiple.R")
source("Final_funcs/gabe_sim_full.R")

start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")

##  We first set the parameters for running in parallel as well as the 
##  population parameter set-up

num.cluster1 = 2
num.runs1 = 1*num.cluster1

k = 1
p = 10
q = 15
Btype = 0
num.obs = 10
n.pair = 2
nperm=100

# cutoff.perc tells where to cutoff for permutation values
cutoff.perc = 0.9

cor.suped = .2
noise.level = 0.1
noise = "clean"
# options are clean, t, sym, and asym  (with sym and asym you need noise.level)
# t uses df=1, we might want a lower df? 1?


B <- build.B(k,p,q,Btype)
#set.seed(47)

run1 <- list()
length(run1)<- num.runs1 #need this val

c1 <- makeCluster(num.cluster1)
registerDoParallel(c1)

##  This loop runs the entire simulation in parallel for each dataset.

perm.results.sim <- foreach(i=1:num.runs1, .combine='cbind') %dopar%{

library(MASS)

simdata = sim.setup.noise.suped(num.obs, B, contamination=noise.level, var.cor=cor.suped, noisetype = noise, Btype=Btype)
sim.output = scca.CVperm(simdata, n.pair, nperm)


# need to check all of this!!
# dimension, quantiles, sorting, etc.

perm.cor.s = sim.output$perm.cor.s
perm.cor.s = apply(t(perm.cor.s), 1, sort, decreasing=T)
perm.s.curve = apply(perm.cor.s, 2, quantile, probs=cutoff.perc)

perm.cor.p = sim.output$perm.cor.p
perm.cor.p = apply(t(perm.cor.p), 1, sort, decreasing=T)
perm.p.curve = apply(perm.cor.p, 2, quantile, probs=cutoff.perc)

res.s = list()
res.s$sp.coef.u = data.frame(matrix(unlist(sim.output$alphas.s), nrow=length(sim.output$alphas.s[[1]]), byrow=F))
res.s$sp.coef.v = data.frame(matrix(unlist(sim.output$betas.s), nrow=length(sim.output$betas.s[[1]]), byrow=F))
res.s$sp.cor = sim.output$cor.test.s

res.p = list()
res.p$sp.coef.u = data.frame(matrix(unlist(sim.output$alphas.p), nrow=length(sim.output$alphas.p[[1]]), byrow=F))
res.p$sp.coef.v = data.frame(matrix(unlist(sim.output$betas.p), nrow=length(sim.output$betas.p[[1]]), byrow=F))
res.p$sp.cor = sim.output$cor.test.p

output.s <- results(res.s, B, n.pair)
output.p <- results(res.p, B, n.pair)


# Is this the output?  Should it be rbind or c?  What do we want to save?
 rbind( interpret.results.curve(output.s, perm.s.curve ),
       interpret.results.curve(output.p, perm.p.curve))


# c(output.jo$cor.test.s, output.jo$perm.cor.s, output.jo$cor.all.s)

}  

fname = paste("Gabecors",Btype,".n",num.obs,".p",p,".q",q,".",noise,".txt",sep="")
write.table(perm.results.sim, file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

end1 <- date()
end1 <- strptime(end1, "%a %b %d %H:%M:%S %Y")
dif1 <- as.numeric(difftime(end1,start,units="mins"))  # how long the first loop takes, in minutes


write.table(cbind(dif1, dif1, fname), file="times.txt", row.names=F, col.names=F, quote=F, sep="\t", append=T)
