
col.names <- c("Num.NonTrival.Pairs","TP.Percent.Tot","TP.Percent.Five","TP.Percent.Cor","Avg.Cor.Tot","Avg.Cor.Five","Num.Vec.Cor",
               "Complete.Groups","Avg.FP.Tot","Avg.FP.Five","Avg.FP.Cor","FN.Tot","FN.Five","FN.Cor")

first <- read.table("Outlier_Sims/ParTrial2.txt")
first.s.bic <- first[,1:13]
first.s.nobic <- first[,14:26]
first.p.bic <- first[,27:39]
first.p.nobic <- first[,40:52]
names(first.s.bic) <- names(first.s.nobic) <- names(first.p.bic) <- names(first.p.nobic) <- col.names


nosvd <- read.table("Outlier_Sims/ParTrial2.NoSVD.txt")
nosvd.s.bic <- nosvd[,1:13]
nosvd.s.nobic <- nosvd[,14:26]
nosvd.p.bic <- nosvd[,27:39]
nosvd.p.nobic <- nosvd[,40:52]
names(nosvd.s.bic) <- names(nosvd.s.nobic) <- names(nosvd.p.bic) <- names(nosvd.p.nobic) <- col.names


noise.01 <- read.table("Outlier_Sims/Noise=.01.2.txt")
noise.01.s.bic <- noise.01[,1:13]
noise.01.s.nobic <- noise.01[,14:26]
noise.01.p.bic <- noise.01[,27:39]
noise.01.p.nobic <- noise.01[,40:52]
names(noise.01.s.bic) <- names(noise.01.s.nobic) <- names(noise.01.p.bic) <- names(noise.01.p.nobic) <- col.names


noise.05 <- read.table("Outlier_Sims/Noise=.05.1.txt")
noise.05.s.bic <- noise.05[,1:13]
noise.05.s.nobic <- noise.05[,14:26]
noise.05.p.bic <- noise.05[,27:39]
noise.05.p.nobic <- noise.05[,40:52]
names(noise.05.s.bic) <- names(noise.05.s.nobic) <- names(noise.05.p.bic) <- names(noise.05.p.nobic) <- col.names


k10 <- read.table("Outlier_Sims/Noise=.01.k=10.txt")
k10.s.bic <- k10[,1:13]
k10.s.nobic <- k10[,14:26]
k10.p.bic <- k10[,27:39]
k10.p.nobic <- k10[,40:52]
names(k10.s.bic) <- names(k10.s.nobic) <- names(k10.p.bic) <- names(k10.p.nobic) <- col.names

eps.5 <- read.table("Outlier_Sims/Eps=.5.txt")
eps.5.s.bic <- eps.5[,1:14]
eps.5.s.nobic <- eps.5[,15:28]
eps.5.p.bic <- eps.5[,29:42]
eps.5.p.nobic <- eps.5[,43:56]
names(eps.5.s.bic) <- names(eps.5.s.nobic) <- names(eps.5.p.bic) <- names(eps.5.p.nobic) <- col.names

