
date()

library(foreach)
library(doParallel)

source("My_Functions/sample.sigma12.function.R")
source("My_Functions/sample.sigma12.function.fullsig.R")
source("My_Functions/bic.function.R")
source("My_Functions/scca.function.R")
source("My_Functions/scca.multiple.R")

num.runs <- 3

G = read.table("CEUsnps.txt")
Y = read.table("CEUexpr.txt")

G = t(G[,-c(1:4)])
Y = t(Y[,-c(1:4)])

p <- dim(G)[2]
q <- dim(Y)[2]
num.obs <- dim(G)[1]

G = matrix(as.numeric(unlist(G)),ncol=p)
Y = matrix(as.numeric(unlist(Y)),ncol=q)


run <- list()
length(run)<- num.runs #need this val

c1 <- makeCluster(10)
registerDoParallel(c1)

perm.results <- foreach(i=1:5, .combine='rbind') %dopar%{

  res.s <- scca.multiple(G, Y, rob = T, bic=T, fullsig=T) #Here bic=T means compute bic as well as nobic
  res.p <- scca.multiple(G, Y, rob = F, bic=T, fullsig=T)

  c(res.s$nobic.cor, res.s$bic.cor, res.p$nobic.cor, res.p$bic.cor)

}  

write.table(perm.results, file="perm.results", row.names=F, quote=F, col.names=F, sep="\t")    

date()