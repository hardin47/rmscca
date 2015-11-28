
# change working directory to desktop rsmcca
library(ggplot2)

n.obs=c(50,100,500,1000)
pq = rbind(c(100,100),c(500,1000))
noise=c("clean","sym","t","asym", "noise")
which.B = c("B1","B2")
full.sim=c()

for(i in c(1:4)){  # num observations
  for (j in 2){ # pq combo
    for(k in c(1,3)){  # noise
      for(l in 2){  # B type

filename = paste("~/Desktop/rsmcca/Final_Sims/secondsim",which.B[l],".n",n.obs[i],".p",pq[j,1],".q",pq[j,2],".",noise[k],".txt", sep="")

sim.new = read.delim(filename, header=FALSE)
names(sim.new) = rep(letters[1:13],2)
reps = dim(sim.new)[1]

# full results

sim.new = cbind(rbind(sim.new[,1:13],sim.new[,14:26]),c(rep("RMSCCA",reps),rep("MSCCA",reps)))


names(sim.new) = c("NNT", "NNT.Perm", "TP.Tot.Perm", "TP.any.Perm", "TP.CG.Perm", 
                   "TP.Rate.Perm", "Ave.Cor.Perm", "CG.Perm", "Avg.FP.Perm", 
                   "Avg.TP.Perm", "FN.Perm", "FDR.Perm", "NCC.Perm","method")
sim.new = cbind(sim.new, "n.obs" = n.obs[i], "p" = pq[j,1], "q" = pq[j,2], "noise"=noise[k], "Btype"=which.B[l])

full.sim = rbind(full.sim, sim.new)

    } }}}


ggplot(full.sim, aes(x=factor(n.obs), y=NNT.Perm, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("number of non-trivial CC")
# + geom_jitter() + coord_trans(y = "sqrt")

pdf(file="TPrateSecond.pdf")
ggplot(full.sim, aes(x=factor(n.obs), y=TP.Rate.Perm, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("true positive rate")
dev.off()

pdf(file="TPCGSecond.pdf")
ggplot(full.sim, aes(x=factor(n.obs), y=TP.CG.Perm/CG.Perm, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("proportion of complete groups")
dev.off()

pdf(file="fdrSecond.pdf")
ggplot(full.sim, aes(x=factor(n.obs), y=FDR.Perm, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("false discovery rate")
dev.off()

pdf(file="fnSecond.pdf")
ggplot(full.sim, aes(x=factor(n.obs), y=FN.Perm/200, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("proportion of false negatives")
dev.off()

pdf(file="nccPairsSecond.pdf")
ggplot(full.sim, aes(x=factor(n.obs), y=NCC.Perm, method,noise)) + 
  geom_boxplot(aes(factor(n.obs), linetype=noise, fill=method)) +
  xlab("number of observations") + ylab("number of significant canonical pairs")
dev.off()



# To calculate power, see which ones are not NA
tapply(!is.na(full.sim$NNT.Perm), list(full.sim$method,full.sim$n.obs,full.sim$p,full.sim$q,full.sim$noise,full.sim$Btype), sum)/100

# TP.Rate.Perm: Spearman is better for noise data, Pearson is better for clean
# TP.Tot.Perm: Spearman is definitely better, maybe even for clean data???
# TP.CG.Perm needs to be divided by CG.Perm 
# NNT.Perm: a little all over the map, possibly Pearson is just generally a lot more canonical pairs
# FDR.Perm: Spearman better, except possibly for clean
# Avg.FP.Perm: Spearman better, except possibly for clean
# TP.any.Perm
# TP.any.Perm / NNT.Perm


# Ave.Cor.Perm: Pearson almost always higher (true generally about Pearson correlation vs. Spearman correlation)
# CG.Perm: ??  Pearson better for clean, Spearman better for noise
# FN.Perm: ?? Spearman might be slightly higher (divide by 48 for B1, 200 for B2)

# Num.Can.Pairs.Perm  -- is that what we should plot for the old figure 2???

