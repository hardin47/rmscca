######################
scca.CVperm <- function(data, n.pair, nperm=100) {
  # Computes cannonical correlations on both the true data and the permuted data
  # Calls sample.sigma12.function and scca.function
  #
  # Is called by fullSimScript.R and nullSimScript.R and dataScript.R
  #
  #
  # Args: 
  #   data - a list containing an n x p matrix called X and an n x q matrix called Y
  #   n.pair - integer, number of canonical pairs to compute
  #
  # Outputs: lam=best.lambdas, cor=best.cor, perm.cor=best.cor.p, cor.all=cor.all, alphas=alphas, betas=betas
  #   lam - matrix of size n.pair by 2 - the cross validated lambdas used at each pair
  #   cor.test - vector of length n.pair - the correlation for each pair computed on the training data
  #   perm.cor - matrix of size 100 by n.pair - the correlation for each permuted pair on the training data
  #   cor.all - vector of length n.pair - correlation on whole data set
  #   alphas - list of loading vectors for X
  #   betas - list of loading vectors for Y
  
  
  n <- nrow(data$X)
  p <- ncol(data$X)
  q <- ncol(data$Y)
  X <- data$X
  Y <- data$Y
  var.root.X <- sqrt(1/diag(var(X)))
  minpq <- min(p,q)  
  Y.old <- Y  
  
  rob.c = c(TRUE, FALSE)
  for (rob in rob.c){
    
    lambda.u <- numeric(n.pair)
    lambda.v <- numeric(n.pair)
    lambda.seq <- seq(0,.5,.1)
    sp.coefs.u.list <- list()
    sp.coefs.v.list <- list()
    sp.coef.u <- c()
    sp.coef.v <- c()
    sp.cor <- c()
    length(sp.coefs.u.list) <- n.pair
    length(sp.coefs.v.list) <- n.pair
    best.cor <- numeric(n.pair)
    best.cor.p <- matrix(0,nrow=nperm, ncol=n.pair)
    n.fold <- 5
    alphas <- list()
    betas <- list()
    cor.all <- c()
    ####Initialization Stuff####
    ####With the permutation####
    for (i.perm in 1:nperm) {
      permute <- sample(1:n)
      Y <- Y[permute,]
      var.root.Y <- sqrt(1/diag(var(Y)))
      k <- list()
      length(k) <- n.pair
      for (i.pair in 1:n.pair) {
        perm <- sample(1:n) #for the cross validation
        n.hold.out <- floor(n / n.fold)
        results <- list()
        for (i.fold in 1:n.fold) {
          results.temp <- c()
          results[[i.fold]] <- matrix(0, nrow=length(lambda.seq)^2, ncol=4)
          test <- perm[ (1 + (i.fold - 1) * n.hold.out) : (i.fold * n.hold.out)]
          test.X <- X[test, ]
          test.Y <- Y[test, ]
          train.X <- X[-test, ]
          train.Y <- Y[-test, ]
          #####Create initial K matrix
          k.train <- sample.sigma12.function(train.X, train.Y, robust = rob) #cross covariance matrix, dimension p x q 
          if (i.pair > 1) {
            for (j.pair in 1:(i.pair-1)) {
              k.train <- k.train - as.vector(t(alphas[[j.pair]]) %*% k.train %*% betas[[j.pair]])  *  alphas[[j.pair]] %*% t(betas[[j.pair]])
            }
          }
          ###initialize cannonical loading vectors
          u.initial.first <- apply(k.train ,1,mean)
          u.initial <- u.initial.first / sqrt(sum(u.initial.first^2))
          v.initial.first <- apply(k.train,2,mean)
          v.initial <- v.initial.first / sqrt(sum(v.initial.first^2))
          
          
          for (l.u in lambda.seq) {
            for (l.v in lambda.seq) {
              uv <- scca.function(k.train, u.initial, v.initial, l.u, l.v) 
              y1 <- test.Y %*% (var.root.Y * uv$v.new)
              x1 <- test.X %*% (var.root.X * uv$u.new)
              y1 <- y1 + (sum(y1 != 0) == 0) * seq(-1, 1, length.out=length(y1))
              x1 <- x1 + (sum(x1 != 0) == 0 ) * y1^2
              if (rob==TRUE){
                cor.temp1 <- cor(x1,y1, method="spearman")
              }else{ cor.temp1 <- cor(x1,y1, method="pearson")}
              results.temp <- rbind(c(i.fold, l.u, l.v, cor.temp1), results.temp)
            }
          }
          results[[i.fold]] <- results.temp  # this is for a fixed permutation
        } #finish i.fold
        
        results.comb <- results[[1]]  #combine results from all folds
        for (i in 2:length(results))
          results.comb <- results.comb + results[[i]]
        results.comb[,4] <- results.comb[,4] / n.fold  #4th column is the correlation
        best <- which.max(results.comb[,4])
        best.lambdas <- results[[1]][best, 2:3]
        best.cor.p[i.perm, i.pair] <- results.comb[best,4]
        
        # Now grab the loadings, so we can update the k matrix, using the whole data
        k[[1]] <- sample.sigma12.function(X, Y, robust = rob)
        if (i.pair > 1) {
          k[[i.pair]] <- k[[i.pair - 1]] - as.vector(t(alphas[[i.pair - 1]]) %*% k[[i.pair - 1]] %*% betas[[i.pair - 1]])  * alphas[[i.pair - 1]] %*% t(betas[[i.pair - 1]])
        }
        u.initial.first <- apply(k[[i.pair]] ,1,mean)
        u.initial <- u.initial.first / sqrt(sum(u.initial.first^2))
        v.initial.first <- apply(k[[i.pair]],2,mean)
        v.initial <- v.initial.first / sqrt(sum(v.initial.first^2))      
        uv.full <- scca.function(k[[i.pair]], u.initial, v.initial, best.lambdas[1], best.lambdas[2])
        alphas[[i.pair]] <- uv.full$u.new
        betas[[i.pair]] <- uv.full$v.new  #these alphas and betas will be used on the folded data 

	if(i.pair==1 & rob==TRUE){lambda1.s = best.lambdas}
	if(i.pair==1 & rob==FALSE){lambda1.p = best.lambdas}
      } # finish i.pair
    } # finish i.perm
    
    
    ####Without the permutation####
    Y <- Y.old
    var.root.Y <- sqrt(1/diag(var(Y)))
    k <- list()
    length(k) <- n.pair
    for (i.pair in 1:n.pair) {
      perm <- sample(1:n) #for the cross validation
      n.hold.out <- floor(n / n.fold)
      results <- list()
      for (i.fold in 1:n.fold) {
        results.temp <- c()
        results[[i.fold]] <- matrix(0, nrow=length(lambda.seq)^2, ncol=4)
        test <- perm[ (1 + (i.fold - 1) * n.hold.out) : (i.fold * n.hold.out)]
        test.X <- X[test, ]
        test.Y <- Y[test, ]
        train.X <- X[-test, ]
        train.Y <- Y[-test, ]
        #####Create initial K matrix
        k.train <- sample.sigma12.function(train.X, train.Y, robust = rob) #cross covariance matrix, dimension p x q 
        if (i.pair > 1) {
          for (j.pair in 1:(i.pair-1)) {
            k.train <- k.train - as.vector(t(alphas[[j.pair]]) %*% k.train %*% betas[[j.pair]])  *  alphas[[j.pair]] %*% t(betas[[j.pair]])
          }
        }
        ###initialize cannonical loading vectors
        u.initial.first <- apply(k.train ,1,mean)
        u.initial <- u.initial.first / sqrt(sum(u.initial.first^2))
        v.initial.first <- apply(k.train,2,mean)
        v.initial <- v.initial.first / sqrt(sum(v.initial.first^2))
        
        for (l.u in lambda.seq) {
          for (l.v in lambda.seq) {
            uv <- scca.function(k.train, u.initial, v.initial, l.u, l.v) 
            y1 <- test.Y %*% (var.root.Y * uv$v.new)
            x1 <- test.X %*% (var.root.X * uv$u.new)
            y1 <- y1 + (sum(y1 != 0) == 0) * seq(-1, 1, length.out=length(y1))
            x1 <- x1 + (sum(y1 != 0) == 0 ) * y1^2
            if (rob==TRUE){
              cor.temp2 <- cor(x1,y1, method="spearman")
            }else{ cor.temp2 <- cor(x1,y1, method="pearson")}
            results.temp <- rbind(c(i.fold, l.u, l.v, cor.temp2), results.temp)
          }
        }
        results[[i.fold]] <- results.temp  # this is for a fixed permutation
      }
      
      results.comb <- results[[1]]  #combine results from all folds
      for (i in 2:length(results))
        results.comb <- results.comb + results[[i]]
      results.comb[,4] <- results.comb[,4] / n.fold  #4th column is the correlation
      best <- which.max(results.comb[,4])
      best.lambdas <- results[[1]][best, 2:3]
      best.cor[i.pair] <- results.comb[best,4]
      
      # Now grab the loadings, so we can update the k matrix, using the whole data
      k[[1]] <- sample.sigma12.function(X, Y, robust = rob)
      if (i.pair > 1) {
        k[[i.pair]] <- k[[i.pair - 1]] - as.vector(t(alphas[[i.pair - 1]]) %*% k[[i.pair - 1]] %*% betas[[i.pair - 1]]) * alphas[[i.pair - 1]] %*% t(betas[[i.pair - 1]])
      }
      u.initial.first <- apply(k[[i.pair]] ,1,mean)
      u.initial <- u.initial.first / sqrt(sum(u.initial.first^2))
      v.initial.first <- apply(k[[i.pair]],2,mean)
      v.initial <- v.initial.first / sqrt(sum(v.initial.first^2))
      current.cor <- -1
      for (i.lam in lambda.seq) {
        for (j.lam in lambda.seq) {
          uv.full <- scca.function(k[[i.pair]], u.initial, v.initial, i.lam, j.lam)
          alphas[[i.pair]] <- uv.full$u.new
          betas[[i.pair]] <- uv.full$v.new  #these alphas and betas will be used on the folded data 
          y1 <- Y %*% (var.root.Y * betas[[i.pair]])
          x1 <- X %*% (var.root.X * alphas[[i.pair]])
          y1 <- y1 + (sum(y1 != 0) == 0) * seq(-1, 1, length.out=length(y1))
          x1 <- x1 + (sum(x1 != 0) == 0 ) * y1^2
          if (rob==TRUE) {
            cor.temp <- cor(x1,y1, method="spearman")
          } else {  
            cor.temp <- cor(x1,y1, method="pearson")
          }
          if (cor.temp > current.cor) {
            current.cor <- cor.temp
            lam.fav <- c(i.lam, j.lam)
          }
        }
      }
      uv.full <- scca.function(k[[i.pair]], u.initial, v.initial, lam.fav[1], lam.fav[2])
      alphas[[i.pair]] <- uv.full$u.new
      betas[[i.pair]] <- uv.full$v.new  #these alphas and betas will be used on the folded data 
      y1 <- Y %*% (var.root.Y * betas[[i.pair]])
      x1 <- X %*% (var.root.X * alphas[[i.pair]])
      y1 <- y1 + (sum(y1 != 0) == 0) * seq(-1, 1, length.out=length(y1))
      x1 <- x1 + (sum(x1 != 0) == 0 ) * y1^2
      if (rob==TRUE){
        cor.all[i.pair] <- cor(x1,y1, method="spearman")
      }else{ cor.all[i.pair] <- cor(x1,y1, method="pearson")}
    } # finish i.pair
    
    if (rob==TRUE){
      best.lam.save.s = lam.fav # lambda on all data (gives very high cor)
      best.cor.s = best.cor     # correlation using test/CV values
      best.cor.p.s = best.cor.p # permutation correlations, to find Q-curve
      cor.all.s = cor.all       # very high cor on *all* data w/o CV
      alphas.s = alphas         # alphas on all data (gives very high cor)
      betas.s = betas           # betas on all data (gives very high cor) 
    }else{
      best.lam.save.p = lam.fav
      best.cor.pears = best.cor
      best.cor.p.p = best.cor.p
      cor.all.p = cor.all
      alphas.p = alphas
      betas.p = betas
    }
    
    
  } # finish robust loop
  
    return(list(lam.s=best.lam.save.s, lambda1.s = lambda1.s, cor.test.s=best.cor.s, perm.cor.s=best.cor.p.s, cor.all.s=cor.all.s, 
                alphas.s=alphas.s, betas.s=betas.s, lam.p=best.lam.save.p, lambda1.p = lambda1.p, cor.test.p=best.cor.pears, 
                perm.cor.p=best.cor.p.p, cor.all.p=cor.all.p, alphas.p=alphas.p, betas.p=betas.p))
    
} # end of function





