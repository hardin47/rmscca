

# running our SCCA analysis on the breastcancer dataset in PMA

date()

library(foreach)
library(doParallel)

# setwd("C:/Users/Jo Hardin/Desktop/rsmcca")
source("Final_funcs/sample_sigma12_function.R")
source("Final_funcs/scca_function.R")
source("Final_funcs/scca_CVperm.R")

start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")


num.cluster1 = 2
num.runs1 = 1*num.cluster1
run1 <- list()
length(run1)<- num.runs1 #need this val

c1 <- makeCluster(num.cluster1)
registerDoParallel(c1)

#p = 2149
#q = 19672
num.obs = 89
n.pair = 10  # should be at least 10
nperm=100
# cutoff.perc tells where to cutoff for permutation values
cutoff.perc = 0.9


# install.packages("PMA", lib = "~/", repos="http://cran.rstudio.com/")
# library(PMA, lib.loc = "~/")
# library(PMA)
# data(breastdata)
# attach(breastdata)

realdata=list()
#realdata$X = t(dna)
#realdata$Y = t(rna)
realdata$dna = read.table("dna.csv", header=F, sep=",")
realdata$rna = read.table("rna.csv", header=F, sep=",")
realdata$genechr = read.table("genechr.csv", header=F, sep=",")
realdata$chrom = read.table("chrom.csv", header=F, sep=",")


comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# running the chromosomes in parallel
results.real <- foreach(i=1:num.runs1, .combine='comb', 
	.init=list(list(), list(), list(), list(),list(), list(), list(), list())) %dopar%{

tempdata = list()
tempdata$X = t(realdata$rna[which(realdata$genechr==i),])
tempdata$Y = t(realdata$dna[which(realdata$chrom==i),])
#tempdata$Y = t(realdata$dna)

real.output = scca.CVperm(tempdata, n.pair, nperm)

# using the permuted correlations to create a curve to determine significance cutoffs
perm.cor.s = real.output$perm.cor.s
perm.s.curve = apply(perm.cor.s, 2, quantile, probs=cutoff.perc)

perm.cor.p = real.output$perm.cor.p
perm.p.curve = apply(perm.cor.p, 2, quantile, probs=cutoff.perc)


# mapping new output to the same form as the previous output
res.s = list()
res.s$sp.coef.u = data.frame(matrix(unlist(real.output$alphas.s), nrow=length(real.output$alphas.s[[1]]), byrow=F))
res.s$sp.coef.v = data.frame(matrix(unlist(real.output$betas.s), nrow=length(real.output$betas.s[[1]]), byrow=F))
res.s$sp.cor = real.output$cor.test.s

res.p = list()
res.p$sp.coef.u = data.frame(matrix(unlist(real.output$alphas.p), nrow=length(real.output$alphas.p[[1]]), byrow=F))
res.p$sp.coef.v = data.frame(matrix(unlist(real.output$betas.p), nrow=length(real.output$betas.p[[1]]), byrow=F))
res.p$sp.cor = real.output$cor.test.p

list(res.s$sp.coef.u, res.s$sp.coef.v, res.s$sp.cor, perm.s.curve, 
	res.p$sp.coef.u, res.p$sp.coef.v, res.p$sp.cor, perm.p.curve)

} 


fname = paste("breastcors.txt",sep="")
write.table(cbind(rbind(results.real[[3]], results.real[[4]]), rbind(results.real[[7]],results.real[[8]])),
		file=fname, row.names=F, quote=F, col.names=F, sep="\t")    

fname = paste("realalphas.txt",sep="")
write.table(cbind(results.real[[1]], results.real[[5]]), file=fname, row.names=F, quote=F, col.names=F, sep="\t")

fname = paste("realbetas.txt",sep="")
write.table(cbind(results.real[[2]], results.real[[6]]), file=fname, row.names=F, quote=F, col.names=F, sep="\t")

#library(fanplot)
#plot(c(1:n.pair), real.output$cor.test.s, ylim=c(0,1), xlab="canonical pair", ylab="correlation")
#fan(real.output$perm.cor.s[-1,], type="percentile", probs=c(0.8,.9,.95,.99))
#points(c(1:n.pair), real.output$cor.test.s, pch=18)
#points(c(1:n.pair), perm.s.curve, pch=25)


#plot(c(1:n.pair), real.output$cor.test.p, ylim=c(0,1), xlab="canonical pair", ylab="correlation")
#fan(real.output$perm.cor.p[-1,], type="percentile", probs=c(0.8,.9,.95,.99))
#points(c(1:n.pair), real.output$cor.test.p, pch=18)
#points(c(1:n.pair), perm.p.curve, pch=25)
