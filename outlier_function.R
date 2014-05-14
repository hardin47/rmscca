source("My_Functions/results.R")
source("My_Functions/interpret_results.R")
source("My_Functions/build_B.R")
source("My_Functions/sim_setup_noise.R")


source("Current_R_Files/Hardin's_Random_Covariance_Function.R")
source("My_Functions/sample.sigma12.function.R")
source("My_Functions/bic.function.R")
source("My_Functions/scca.function.R")
source("My_Functions/scca.multiple.R")
source("My_Functions/select.parameters.multiple.R")

num.runs <- 35

k <- .95
p <- 750
q <- 30
num.obs <- 400

cor.val.cutoff <- .4

#Couldn't parameterize the B structure. It's hardcoded in the function
B <- build.B(k,p,q)

run <- list()
length(run)<- num.runs #need this val

for(i in 1:num.runs){
  setup <- sim.setup(num.obs,B)
  G <- setup$X
  Y <- setup$Y
  
  start <- date()
  start <- strptime(start, "%a %b %d %H:%M:%S %Y")
  res.s <- scca.multiple(G,Y, rob = T, bic=T) #Here bic=T means compute bic as well as nobic
  res.p <- scca.multiple(G,Y, rob = F, bic=T)
  end <- date()
  end <- strptime(end, "%a %b %d %H:%M:%S %Y")
  difftime(end, start, units = "mins")
  
  output.s.bic <- results(res.s, B, bic=T)
  output.s.nobic <- results(res.s, B, bic=F)
  
  output.p.bic <- results(res.p, B, bic = T)
  output.p.nobic <- results(res.p, B, bic = F)
  
  big.s.bic <- big.s.nobic <- big.p.bic <- big.p.nobic <- matrix(nrow = num.runs, ncol = 13)
  
  big.s.bic[i,] <-  interpret.results(output.s.bic,cor.val.cutoff)
  big.s.nobic[i,] <-interpret.results(output.s.nobic,cor.val.cutoff)
  big.p.bic[i,] <- interpret.results(output.p.bic,cor.val.cutoff)
  big.p.nobic[i,] <- interpret.results(output.p.nobic, cor.val.cutoff)
  
  big.s.bic <- as.data.frame(big.s.bic)
  big.s.nobic <- as.data.frame(big.s.nobic)
  big.p.bic <- as.data.frame(big.p.bic)
  big.p.nobic <- as.data.frame(big.p.nobic)
 
  
  col.names <- c("Num.NonTrival.Pairs","TP.Percent.Tot","TP.Percent.Five","TP.Percent.Cor","Avg.Cor.Five","Num.Vec.Cor",
             "Complete.Groups","Avg.FP.Tot","Avg.FP.Five","Avg.FP.Cor","FN.Tot","FN.Five","FN.Cor")
  
  names(big.s.bic)<-names(big.s.nobic)<-names(big.p.bic)<-names(big.p.nobic) <- col.names
  
  big <- matrix(nrow=num.runs,ncol=(12*4))
  View(c(interpret.results(output.s.bic,cor.val.cutoff),
    interpret.results(output.s.nobic,cor.val.cutoff),
    interpret.results(output.p.bic,cor.val.cutoff),
    interpret.results(output.p.nobic, cor.val.cutoff)))
  
    
}
