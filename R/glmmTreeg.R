######=====================================================================================
###### functions for continous outcome
######=====================================================================================

CCov_treekernel <- function(D, alpha){
  exp(-(D^2)*alpha)
}


glmmTreeg <-function(y,Z,X=NULL,D,lambda1,lambda2,norm.kernel=F){
  
  nlambda1.0=nlambda1=length(lambda1)
  nlambda2.0=nlambda2=length(lambda2)
  lambda1.0=lambda1
  lambda2.0=lambda2
  
  n=length(y)
  p.z=ncol(Z)

  if(is.null(X)){
    XX=as.matrix(rep(1,n))
  }else{
    XX=cbind(1,X)
  }
  
  mind=matrix(1,nlambda1,nlambda2)
  dimnames(mind)=list(lambda1=lambda1,lambda2=lambda2)
  obj=list()
  
  for(i in 1:nlambda1){
    ZZ=Z^lambda1[i]
    ZZ[Z==0]=0
    for(j in 1:nlambda2){
      K=CCov_treekernel(D=D, alpha=lambda2[j])
      K=as.matrix(nearPD(K)$mat)
      if(norm.kernel) K=cov2cor(K) 
      er1=try(obj[[paste(lambda1[i])]][[paste(lambda2[j])]]<-mixed.solve(y=y, Z=ZZ, K=K, X=XX))	
      if(er1[1]=="Error in mixed.solve(y = y, Z = ZZ, K = K, X = Xx, method = Method, bounds = bounds,  : \n  K not positive semi-definite.\n"){
        K=nearPD(K)$mat
        er2=try(obj[[paste(lambda1[i])]][[paste(lambda2[j])]]<-mixed.solve(y=y, Z=ZZ, K=K, X=XX))
        if(er2[1]=="Error in diag(K) : no method for coercing this S4 class to a vector\n"){
          mind[i,j]=0
        }
      }	
    }
  }
  
  
  lambda1=lambda1[apply(mind,1,sum)!=0]
  lambda2=lambda2[apply(mind,2,sum)!=0]
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  mind=as.matrix(mind[apply(mind,1,sum)!=0,apply(mind,2,sum)!=0,drop=FALSE])

  if(nlambda1>0){
    
    us=betas=list()
    Ls=matrix(NA,nlambda1,nlambda2)
    dimnames(Ls)=list(lambda1=lambda1,lambda2=lambda2)
    for(i in 1:nlambda1){
      beta=u=NULL
      for(j in 1:nlambda2){
        if(mind[i,j]==1){
          beta=cbind(beta,as.matrix(obj[[paste(lambda1[i])]][[paste(lambda2[j])]]$beta))
          u=cbind(u,as.matrix(obj[[paste(lambda1[i])]][[paste(lambda2[j])]]$u))
          Ls[i,j]=obj[[paste(lambda1[i])]][[paste(lambda2[j])]]$LL
        }
      }
      if(is.null(X)){
        dimnames(beta)=list('Intercept', round(lambda2.0[mind[i,]!=0],digits = 5))
        betas[[paste(lambda1[i])]]=beta
      }else{
        beta0=matrix(0, nrow = ncol(XX), ncol = length(lambda2.0[mind[i,]!=0]))
        beta0[1,]=beta[1,]
        beta0[2:ncol(XX),]=beta[-1,]
        varnames=colnames(XX)
        if (is.null(colnames(X))) 
          varnames=paste("V", 1:ncol(X), sep = "")
        varnames=c("Intercept", varnames)
        dimnames(beta0)=list(varnames, round(lambda2.0[mind[i,]!=0],digits = 5))
        betas[[paste(lambda1[i])]]=beta0
      }
      
      varnames=colnames(Z)
      if(is.null(colnames(Z))) 
        varnamesppaste("Z", 1:ncol(Z), sep = "")
      dimnames(u)=list(varnames, round(lambda2.0[mind[i,]!=0],digits = 5))
      us[[paste(lambda1[i])]]=u
    }
    
    names(betas)=names(us)=lambda1
  }else{
      warning("No solution for current parameter setting!")
      return(NULL)
  }
  
  list(betas=betas,us=us,Ls=Ls,lambda1=lambda1,lambda2=lambda2,nlambda1=nlambda1,nlambda2=nlambda2,mind=mind,lambda1.0=lambda1.0,lambda2.0=lambda2.0)
  
}


predict.glmmTreeg<-function(obj,X,Z){
    
    lambda1=obj$lambda1
    lambda2=obj$lambda2
    nlambda1=length(lambda1)
    nlambda2=length(lambda2)
    etas=list()
    for(i in 1:nlambda1){
        ZZ=Z^lambda1[i]
        ZZ[Z==0]=0
        etas[[paste(lambda1[i])]]=as.matrix(ZZ %*% obj$us[[paste(lambda1[i])]])
        if(is.null(X)){
          etas[[paste(lambda1[i])]]=sweep(etas[[paste(lambda1[i])]],2,obj$betas[[paste(lambda1[i])]],"+")
        }else{
          etas[[paste(lambda1[i])]]=sweep(etas[[paste(lambda1[i])]]+X %*% obj$betas[[paste(lambda1[i])]][-1,,drop=FALSE],2,obj$betas[[paste(lambda1[i])]][1,,drop=FALSE],"+")
        }
    }
    
    etas
    
  }
  
  
  getmin2 <- function(lambda,cvm,cvsd){
    cvmin <- min(cvm,na.rm=TRUE)
    idmin <- cvm<=cvmin
    lambda.min <- max(lambda[idmin],na.rm=TRUE)
    list(lambda.min=lambda.min,cvmin=cvmin)
  }
  

  setcritcv <- function(crit,cv.ind,lambda1,lambda2,nfolds){

    nlambda1=length(lambda1)
    nlambda2=length(lambda2)
    outmat <-  array(NA,c(nfolds,nlambda1,nlambda2))
    good  <- array(0,c(nfolds,nlambda1,nlambda2))
    dimnames(outmat)=dimnames(good)=list(nfolds=1:nfolds,lambda1=lambda1,lambda2=lambda2)

    for(ifold in seq(nfolds)){
      whichi <- cv.ind==ifold
      for(i in seq(nlambda1)){  # different #lambda2 for different fold
        outmat[ifold,i,]=apply(crit[whichi,i,,drop=FALSE],c(2,3),mean,na.rm=TRUE)
        good[ifold,i,][!is.na(outmat[ifold,i,])]=1
      }
    }
    
    Nmat <- matrix(0, nlambda1,nlambda2)
    for(i in seq(nlambda1)){
        for(j in seq(nlambda2)){
          Nmat[i,j] <- sum(good[,i,j])
        }
    }
    
    list(crit=outmat,Nmat=Nmat)
  }
  
  

  setcrit <- function(obj,X,Z,y,crit,whichi){
    pred <- predict.glmmTreeg(obj, X,Z)
    for(i in seq(obj$nlambda1)) crit[whichi,i,seq(ncol(pred[[i]]))] <- ( sweep( pred[[i]],1,y,FUN="-") )^2
    crit
  }
  
  
  
  
  predict.cv.glmmTreeg<-function(obj.cv,X,Z){
    predict.glmmTreeg(obj.cv$obj.min,X,Z)[[1]]
  }
  
  
  


cv.glmmTreeg<-function(y,Z,X=NULL,D,lambda1,lambda2, nfolds=5,trace=TRUE,norm.kernel=F){
  
  obj=glmmTreeg(y,Z,X,D,lambda1,lambda2,norm.kernel=norm.kernel)
  if (is.null(obj)) 
    stop("No solution for current parameter setting!")
  
  n=length(y)
  p.z=ncol(Z)
  
  
  lambda1=obj$lambda1
  lambda2=obj$lambda2
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  mind=obj$mind
  
  
  crit <- array(NA, c(n,nlambda1, nlambda2), dimnames=list(y=1:n, paste(lambda1),paste(lambda2)))
  cv.ind <- ceiling(sample(1:n)/n*nfolds)
  
  cv.args=list()
  cv.args$y=y
  cv.args$Z=Z
  cv.args$X=X
  cv.args$D=D
  cv.args$lambda1=lambda1
  cv.args$lambda2=lambda2
  
  
  for (ifold in 1:nfolds) {
    if(trace) cat("Starting CV fold #",ifold,sep="","\n")
    whichi <- cv.ind==ifold
    cv.args$X <- X[!whichi, , drop=FALSE]
    cv.args$y <- y[!whichi]
    cv.args$Z <- Z[!whichi, , drop=FALSE]
    
    obj.i <- do.call(glmmTreeg,cv.args)
    
    if(is.null(obj.i)) stop("No solution within nfold, please reset your folds again!")
    X2 <- X[whichi, , drop=FALSE]
    Z2 <- Z[whichi, , drop=FALSE]
    y2 <- y[whichi]
    crit=setcrit(obj.i,X2,Z2,y2,crit,whichi)
  }
  
  
  critcv=setcritcv (crit,cv.ind,lambda1,lambda2,nfolds)
  
  
  ## Return
  cvm <- cvsd <- list()
  cvmat <- data.frame(matrix(rep(0,3*nlambda1),nrow=nlambda1))
  colnames(cvmat) <- c("lambda1","lambda2.min","cvmin")
  
  for(i in 1:nlambda1) {
    ind <- critcv$Nmat[i,] >=nfolds
    if(all(ind==FALSE)) {
      warning("No criteria value for all folds for some lambda2, may need to adjust tuning parameters.")
      next;		
    }
    cvm[[i]] <- apply(critcv$crit[,i,ind,drop=FALSE],c(2,3),mean,na.rm=TRUE)
    cvsd[[i]] <- sqrt( apply(scale(critcv$crit[,i,ind],cvm[[i]],FALSE)^2,2,mean,na.rm=TRUE)/critcv$Nmat[i,ind])
    lmin <- getmin2(as.numeric(colnames(cvm[[i]])),cvm[[i]],cvsd[[i]])
    cvmat[i,] <- c(lambda1[i],lmin$lambda.min, lmin$cvmin)
  } 
  
  cvm <- cvm[!sapply(cvm, is.null)]
  cvsd <- cvsd[!sapply(cvsd, is.null)]
  cvmat <- cvmat[!apply(cvmat, 1, sum) == 0, ]
  id.min <- which.min(cvmat$cvmin)
  cvmin <- min(cvmat$cvmin)
  lambda1.min <- cvmat$lambda1[id.min]
  lambda2.min <- cvmat$lambda2.min[id.min]
  cv.args$X <- X
  cv.args$y <- y
  cv.args$Z <- Z
  cv.args$lambda1 <- lambda1.min
  cv.args$lambda2 <- lambda2.min
  obj.min <- do.call(glmmTreeg, cv.args)
  
  
  obj.cv <- list(obj.min=obj.min,cvmat=cvmat,cvmin=cvmin,cvm=cvm,cvsd=cvsd,lambda1.min=lambda1.min,lambda2.min=lambda2.min,
                 beta.min=obj.min$betas[[1]],u.min=obj.min$us[[1]])
  obj.cv
  
  
}









