

# running our SCCA analysis on the breastcancer dataset in PMA
#  ONLY chromosome 2

date()

#####  get functions
# setwd("C:/Users/Jo Hardin/Desktop/rsmcca")
source("Final_funcs/sample_sigma12_function.R")
source("Final_funcs/scca_function.R")
source("Final_funcs/scca_CVperm.R")



#### set parameter values
num.obs = 89
n.pair = 10  # should be at least 10
nperm=100
# cutoff.perc tells where to cutoff for permutation values
cutoff.perc = 0.9


#### input data
# library(PMA)
# data(breastdata)
# attach(breastdata)
realdata$dna = read.table("dna.csv", header=F, sep=",")
realdata$rna = read.table("rna.csv", header=F, sep=",")
realdata$genechr = read.table("genechr.csv", header=F, sep=",")
realdata$chrom = read.table("chrom.csv", header=F, sep=",")

tempdata = list()
tempdata$X = t(realdata$rna[which(realdata$genechr==2),])
tempdata$Y = t(realdata$dna[which(realdata$chrom==2),])


#### Permutations
# running the permutation piece to get the appropriate cutoff values
# the permutation (including CV) must be run on the actual data so that 
# the cutoffs are appropriate for the data at hand
 
real.output = scca.CVperm(tempdata, n.pair, nperm)

# using the permuted correlations to create a curve to determine significance cutoffs
perm.cor.s = real.output$perm.cor.s
perm.cor.s = apply(t(perm.cor.s), 1, sort, decreasing=T)
perm.s.curve = apply(perm.cor.s, 2, quantile, probs=cutoff.perc)

perm.cor.p = real.output$perm.cor.p
perm.cor.p = apply(t(perm.cor.p), 1, sort, decreasing=T)
perm.s.curve = apply(perm.cor.p, 2, quantile, probs=cutoff.perc)

# mapping new output to the same form as the previous output
res.s = list()
res.s$sp.coef.u = data.frame(matrix(unlist(real.output$alphas.s), nrow=length(real.output$alphas.s[[1]]), byrow=F))
res.s$sp.coef.v = data.frame(matrix(unlist(real.output$betas.s), nrow=length(real.output$betas.s[[1]]), byrow=F))
res.s$sp.cor = real.output$cor.test.s

res.p = list()
res.p$sp.coef.u = data.frame(matrix(unlist(real.output$alphas.p), nrow=length(real.output$alphas.p[[1]]), byrow=F))
res.p$sp.coef.v = data.frame(matrix(unlist(real.output$betas.p), nrow=length(real.output$betas.p[[1]]), byrow=F))
res.p$sp.cor = real.output$cor.test.p

results.real = list(res.s$sp.coef.u, res.s$sp.coef.v, res.s$sp.cor, perm.s.curve, res.p$sp.coef.u, res.p$sp.coef.v, res.p$sp.cor, perm.p.curve)

 


fname = paste("breastcors2.txt",sep="")
write.table(cbind(rbind(results.real[[3]], results.real[[4]]), rbind(results.real[[7]],results.real[[8]])),
		file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

fname = paste("realalphas2.txt",sep="")
write.table(cbind(results.real[[1]], results.real[[5]]), file=fname, row.names=F, quote=F, col.names=F, sep="\t")

fname = paste("realbetas2.txt",sep="")
write.table(cbind(results.real[[2]], results.real[[6]]), file=fname, row.names=F, quote=F, col.names=F, sep="\t")

date()




