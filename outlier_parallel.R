library(foreach)
library(doParallel)

source("My_Functions2/results.R")
source("My_Functions2/interpret_results.R")
source("My_Functions2/build_B.R")
source("My_Functions2/sim_setup_noise.R")
source("My_Functions2/sim_setup.R")
source("My_Functions2/sim_setup_suped.R")
source("My_Functions2/sim_setup_noise_suped.R")


source("Current_R_Files/Hardin's_Random_Covariance_Function.R")
source("My_Functions2/randCov_suped.R")
source("My_Functions2/sample.sigma12.function.R")
source("My_Functions2/bic.function.R")
source("My_Functions2/scca.function.R")
source("My_Functions2/scca.multiple.R")
source("My_Functions2/select.parameters.multiple.R")


start <- date()
start <- strptime(start, "%a %b %d %H:%M:%S %Y")

num.cluster <- 30
num.runs <- num.cluster

k <- .95
p <- 750
q <- 30
num.obs <- 400

cor.val.cutoff <- .4
cor.suped <- .3
noise.level <- .05

#Couldn't parameterize the B structure. It's hardcoded in the function
B <- build.B(k,p,q)

c1 <- makeCluster(num.cluster)
registerDoParallel(c1)

sim.results <- foreach(i=1:num.runs, .combine='rbind') %dopar%{

  setup <- sim.setup.noise.suped(num.obs,B,contamination = noise.level, var.cor  = cor.suped)
  G <- setup$X
  Y <- setup$Y
  res.s <- scca.multiple(G,Y, rob = T, bic=T) #Here bic=T means compute bic as well as nobic
  res.p <- scca.multiple(G,Y, rob = F, bic=T)
  
  
  output.s.bic <- results(res.s, B, bic=T)
  output.s.nobic <- results(res.s, B, bic=F)
  
  output.p.bic <- results(res.p, B, bic = T)
  output.p.nobic <- results(res.p, B, bic = F)
  
#   big.s.bic <- big.s.nobic <- big.p.bic <- big.p.nobic <- matrix(nrow = num.runs, ncol = 12)
#   
#   big.s.bic[i,] <-  interpret.results(output.s.bic,cor.val.cutoff)
#   big.s.nobic[i,] <-interpret.results(output.s.nobic,cor.val.cutoff)
#   big.p.bic[i,] <- interpret.results(output.p.bic,cor.val.cutoff)
#   big.p.nobic[i,] <- interpret.results(output.p.nobic, cor.val.cutoff)
  
 c(interpret.results(output.s.bic,cor.val.cutoff),
    interpret.results(output.s.nobic,cor.val.cutoff),
    interpret.results(output.p.bic,cor.val.cutoff),
    interpret.results(output.p.nobic, cor.val.cutoff))
}  
    

write.table(sim.results, file = "noise.05.suped.B.Eps=.5.txt",quote=F,col.names=F,row.names=F)
  end <- date()
  end <- strptime(end, "%a %b %d %H:%M:%S %Y")
  dif <- as.numeric(difftime(end, start, units = "mins"))

write.table(dif, file = "TimeDif.noise.05.suped.B.Eps=.5.txt")


