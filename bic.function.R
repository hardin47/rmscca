bic.function<-function(uv,X,Y,p,q,n.sample,robust=T){
  
  uv.bound<-as.matrix(cbind(t(uv$u.new),t(uv$v.new)))
  uv.bound.loop<-uv.bound
  d<-sum(abs(uv.bound)>0)
  index<-array(dim=(d-1))
  r<-array(dim=(d-1))
  BIC<-array(dim=(d-1))

  cov.x = cov(X)
  cov.y = cov(Y)
  diag(cov.x) = ifelse(diag(cov.x) == 0 , 1, diag(cov.x))
  diag(cov.y) = ifelse(diag(cov.y) == 0, 1, diag(cov.y))
  
  if(robust==T){
    r[1]<-cor(X%*%diag(sqrt(1/diag(cov.x)))%*%uv.bound.loop[1:p],Y%*%diag(sqrt(1/diag(cov.y)))%*%uv.bound.loop[(p+1):(p+q)], method="spearman")
  }else{
    r[1]<-cor(X%*%diag(sqrt(1/diag(cov.x)))%*%uv.bound.loop[1:p],Y%*%diag(sqrt(1/diag(cov.y)))%*%uv.bound.loop[(p+1):(p+q)], method="pearson")
  }   
  BIC[1]<-n.sample*log(1-r[1]^2)+(d)*log(n.sample)
  index[1]<-0
  
  if(d>2){
    for (j in 2:(d-1)){
      if (sum(uv.bound.loop[1:p] != 0) == 1){
        uv.temp<-uv.bound.loop[(p+1):(p+q)]
        uv.min = min(uv.temp[uv.temp!=0])
        i<-which(uv.bound.loop == uv.min)
        index[j]<-i[1]
        uv.bound.loop[i] <- 0
      }
      else if (sum(uv.bound.loop[(p+1):(p+q)] != 0) == 1){
        uv.temp<-uv.bound.loop[1:p]
        uv.min = min(uv.temp[uv.temp!=0])
        i<-which(uv.bound.loop == uv.min)
        index[j]<-i[1]
        uv.bound.loop[i] <- 0
      }else{
        uv.min = min(uv.bound.loop[uv.bound.loop!=0])
        i<-which(uv.bound.loop == uv.min)
        index[j]<-i[1]
        uv.bound.loop[i] <- 0
      }
      if(robust==T){      
        r[j]<-cor(X%*%diag(sqrt(1/diag(cov.x)))%*%uv.bound.loop[1:p],Y%*%diag(sqrt(1/diag(cov.y)))%*%uv.bound.loop[(p+1):(p+q)], method="spearman")
      }else{
        r[j]<-cor(X%*%diag(sqrt(1/diag(cov.x)))%*%uv.bound.loop[1:p],Y%*%diag(sqrt(1/diag(cov.y)))%*%uv.bound.loop[(p+1):(p+q)], method="pearson")
      }
      BIC[j]<-n.sample*log(1-r[j]^2)+(d-(j-1))*log(n.sample)
    }
  }
  
  k<-which.min(BIC)
  if(BIC[k]>0){
    uv.bound <- rep(0,length(uv.bound))
  }else{
    uv.bound[index[1:k]]<-0
  }
  
  #At this point, we may want to run CCA one more time on the nonzero variables to adjust the coefficients and improve the correlation.
  return(list(coefs.u = uv.bound[1:p],coefs.v = uv.bound[(p+1):(p+q)]))
}
