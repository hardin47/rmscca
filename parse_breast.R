
# using the data from the PMA package on breast cancer, we use MSCCA to look
# at the linear realationships for each chromosome.  

####################################################################
## Entering the data from the simulation
####################################################################

setwd("C:/Users/Jo Hardin/Desktop/rsmcca")
breastcor = read.table("Final_Sims/breastcors.txt", header=F, sep="\t")
breastalpha = read.table("Final_Sims/breastalphas.txt", header=F, sep="\t")
breastbeta = read.table("Final_Sims/breastbetas.txt", header=F, sep="\t")
library(PMA)
data(breastdata)
attach(breastdata)

num.runs1 = 23  # number of parallel runs, here number of chrom
n.pair=10

breastcor.s = breastcor[,1:(num.runs1*n.pair)]
breastcor.p = breastcor[,(num.runs1*n.pair+1):(2*num.runs1*n.pair)]

bcor.obs.s = array(dim=c(num.runs1, n.pair))
bcor.q.s = array(dim=c(num.runs1, n.pair))
bcor.obs.p = array(dim=c(num.runs1, n.pair))
bcor.q.p = array(dim=c(num.runs1, n.pair))

for(i in 1:num.runs1){
bcor.obs.s[i,] = unlist(breastcor.s[1,((i-1)*n.pair + 1):(i*n.pair)])
bcor.q.s[i,] = unlist(breastcor.s[2,((i-1)*n.pair + 1):(i*n.pair)])
bcor.obs.p[i,] = unlist(breastcor.p[1,((i-1)*n.pair + 1):(i*n.pair)])
bcor.q.p[i,] = unlist(breastcor.p[2,((i-1)*n.pair + 1):(i*n.pair)])
}

alpha.s = list()
alpha.p = list()
beta.s = list()
beta.p = list()

for(i in 1){
  alpha.s[[i]] = 
    matrix(breastalpha[1:(cumsum(table(genechr)[1:i])[i]*n.pair),1], ncol=n.pair) 
  alpha.p[[i]] = 
    matrix(breastalpha[1:(cumsum(table(genechr)[1:i])[i]*n.pair),2], ncol=n.pair) 
  
  beta.s[[i]] = 
    matrix(breastbeta[1:(cumsum(table(chrom)[1:i])[i]*n.pair),1], ncol=n.pair) 
  beta.p[[i]] = 
    matrix(breastbeta[1:(cumsum(table(chrom)[1:i])[i]*n.pair),2], ncol=n.pair) 
}
for(i in 2:num.runs1){
  alpha.s[[i]] = 
    matrix(breastalpha[(cumsum(table(genechr)[1:(i-1)])[i-1]*n.pair+1):(cumsum(table(genechr)[1:i])[i]*n.pair),1], ncol=n.pair) 
  alpha.p[[i]] = 
    matrix(breastalpha[(cumsum(table(genechr)[1:(i-1)])[i-1]*n.pair+1):(cumsum(table(genechr)[1:i])[i]*n.pair),2], ncol=n.pair) 
  
  beta.s[[i]] = 
    matrix(breastbeta[(cumsum(table(chrom)[1:(i-1)])[i-1]*n.pair+1):(cumsum(table(chrom)[1:i])[i]*n.pair),1], ncol=n.pair) 
  beta.p[[i]] = 
    matrix(breastbeta[(cumsum(table(chrom)[1:(i-1)])[i-1]*n.pair+1):(cumsum(table(chrom)[1:i])[i]*n.pair),2], ncol=n.pair)   
}


####################################################################
## Finding the significant canonical pairs
####################################################################

numsig.s = c()
numsig.p = c()

for(i in 1:num.runs1){
numsig.s[i] = min(which(bcor.q.s[i,] > bcor.obs.s[i,])) - 1
numsig.p[i] = min(which(bcor.q.p[i,] > bcor.obs.p[i,])) - 1
}

temp.q = bcor.q.s[,1]
temp.obs = bcor.obs.s[,1]

temp = cbind(numsig.s, numsig.p, table(chrom), table(genechr)[1:23])
write.table(temp, "numbersig.txt", sep="\t")
for(i in 1:num.runs1){
write.table(alpha.s[[i]], paste("alpha",i,"s.txt", sep=""), sep="\t")
write.table(beta.s[[i]], paste("beta",i,"s.txt", sep=""), sep="\t")
}



i=1
plot(c(1:10), bcor.obs.p[i,], pch=18, ylim=c(0,1))
points(c(1:10), bcor.q.p[i,], pch=25)
print(i);i=i+1

for (i in 1:23){
  print(c(i, apply(alpha.s[[i]]!=0,2,sum)))
  print(c(i, apply(beta.s[[i]]!=0,2,sum)))
}

######################################################
## Spearman analysis, Chrom 8
######################################################

alpha8 = alpha.s[[8]]  # 673 x 10
beta8 = beta.s[[8]]    # 138 x 10

alpha8.gns = genenames[genechr==8 & !is.na(genechr)]

######################################################
## Real data cors
######################################################

## by chrom

q.firstpair = bcor.q.s[,1]
test.firstpair = bcor.obs.s[,1]

train.firstpair = c()
for (i in 1:num.runs1){
  train.firstpair[i] = cor(t(t(alpha.s[[i]][,1]) %*% rna[genechr==i & !is.na(genechr),]),
                           t(t(beta.s[[i]][,1]) %*% dna[chrom==i & !is.na(chrom),]), method="spearman")
}

plot(q.firstpair, type="b", pch=18, lty=1, xlab="Chromosome", ylab="correlation", ylim=c(0,1.2))
points(test.firstpair, type="b", pch=18, lty=1, col="red")
points(train.firstpair, type="b", pch=18, lty=1, col="purple")
legend(x=18,y=1.3,bty="n",col=c("purple", "red", "black"), pch=18, legend=c("training", "test", "permuted"))


# by pair
# chrom 2

q.chrom2 = bcor.q.s[2,]
test.chrom2 = bcor.obs.s[2,]

train.chrom2 = c()
for (i in 1:n.pair){
  train.chrom2[i] = cor(t(t(alpha.s[[2]][,i]) %*% rna[genechr==2 & !is.na(genechr),]),
                           t(t(beta.s[[2]][,i]) %*% dna[chrom==2 & !is.na(chrom),]), method="spearman")
}

plot(q.chrom2, type="b", pch=18, lty=1, xlab="Canonical Pair", ylab="correlation", 
      main = "Chromosome 2", ylim=c(0,1.2))
points(test.chrom2, type="b", pch=18, lty=2, col="red")
points(train.chrom2, type="b", pch=18, lty=3, col="purple")
legend(x=8,y=1.3,bty="n",col=c("purple", "red", "black"), pch=18, pt.cex=2, lty=3:1, legend=c("training", "test", "permuted"))
