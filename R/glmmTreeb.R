######=====================================================================================
###### functions for binary outcome
######=====================================================================================


logit <- function (x) {
  log(x / (1 - x))
}
inv.logit <- function (x) {
  exp(x) / (1 + exp(x))
}	



cor_tree <- function(value, paras, form = ~1, fixed = TRUE) {
  
  # Place holder - check the validity of the parameter
  object <- value
  attr(object, "formula") <- form
  attr(object, "fixed") <- fixed
  attr(object, "paras") <- paras
  class(object) <- c("cor_tree", "corStruct")
  object
}

Initialize.cor_tree <- function (object, data, ...) {
  
  # Place holder - check the validity of the parameter
  form <- formula(object)
  if (!is.null(getGroupsFormula(form))) {
    attr(object, "groups") <- getGroups(object, form, data = data)
    attr(object, "Dim") <- Dim(object, attr(object, "groups"))
  } else {
    attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
  }
  attr(object, "covariate") <- getCovariate(object, data = data)
  object
}

# Actual computation of the correlation matrix
corMatrix.cor_tree <- function (object, covariate = getCovariate(object), ...) {
  
  # Add an identity matrix 
  paras <- attr(object, "paras")
  hyper <- as.vector(object)
  
  V <- paras[['V']]
  W <- paras[['W']]
  
  if (hyper>30)
  {hyper <- 30}
  diag(V) <- exp(hyper) * diag(V) + W	
  V <- nearPD(V)$mat
  V
}

# Extract the coefficient
coef.cor_tree <- function (object,  ...) {
  
  paras <- attr(object, "paras")
  coefs <- as.vector(object)	
  coefs <- exp(coefs)	
  names(coefs) <- paste("Hyper", 1:length(coefs), sep="")
  coefs
}	


# The code is modified on the glmmPQL from MASS package
glmmTree <- function (model, family, data, weights, 
                      lambda1, lambda2, lambda3=2, D, Z, vc.init=0, norm.kernel=TRUE,
                      control, niter = 10, verbose = TRUE, ...) {

  if (!require("nlme")) 
    stop("package 'nlme' is essential")
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  m <- mcall <- Call <- match.call()
  nm <- names(m)[-1L]
  keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
  for (i in nm[!keep]) m[[i]] <- NULL
  allvars <- all.vars(model) 
  
  Terms <- if (missing(data)) 
    terms(model)
  else terms(model, data = data)
  off <- attr(Terms, "offset")
  if (length(off <- attr(Terms, "offset"))) 
    allvars <- c(allvars, as.character(attr(Terms, "variables"))[off + 1])
  

  #######################################
  # Calculate V
  ZZ <- Z ^ lambda1
  ZZ[Z==0] <- 0
  ZZ <- scale(ZZ, scale=F)
  Zm <- attr(ZZ, 'scaled:center')
  V.tree <- exp(-D^lambda3 * lambda2)
  V <-  ZZ %*%  V.tree %*% t(ZZ) 
  
  if (norm.kernel) {
    Zs <- 1/sqrt(diag(V))
    V <- cov2cor(V)
  }
  
  #######################################
  
  Call$model <- eval(model)
  m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
  environment(m$formula) <- environment(model)
  m$drop.unused.levels <- TRUE
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(m)
  off <- model.offset(mf)
  if (is.null(off)) 
    off <- 0
  wts <- model.weights(mf)
  if (is.null(wts)) 
    wts <- rep(1, nrow(mf))
  mf$wts <- wts
  fit0 <- glm(formula = model, family = family, data = mf, 
              weights = wts, ...)
  w <- fit0$prior.weights
  eta <- fit0$linear.predictors
  yy <- eta + fit0$residuals - off
  wz <- fit0$weights
  fam <- family
  nm <- names(mcall)[-1L]
  keep <- is.element(nm, c("model", "data", "subset", 
                           "na.action", "control"))
  for (i in nm[!keep]) mcall[[i]] <- NULL
  model[[2L]] <- quote(yy)
  mcall[["model"]] <- model
  mcall[[1L]] <- quote(nlme::gls)
  mcall$method <- "ML"
  
  mf$yy <- yy
  mf$wts <- wz
  mf$invwt <- 1/wz
  #######################################
  # Specify correlation(covariance) matrix
  mcall$correlation <- cor_tree(value = vc.init, paras = list(V=V, W=mf$invwt))
  #######################################
  mcall$data <- mf
  
  for (i in seq_len(niter)) {
    
    # updated
    if (verbose) 
      err <- try(fit <- eval(mcall))
    if (class(err) == 'try-error') {
      fit <- list(converge = FALSE)
      return(fit)
    }
    
    #######################################
    # Modified based on PQL
    # hest: Microbiome effect
    hyper <- coef(fit$modelStruct$corStruct)
    V2 <- V1 <- V 
    diag(V2) <- hyper * diag(V2) +  mf$invwt
    hest <- as.vector(V1 %*% solve(V2) %*% residuals(fit, type='response'))
    ########################################		
    etaold <- eta
    eta <- fitted(fit) + off
    if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) 
      break
    
    #######################################
    # Adding the random effects
    etab <- eta + hest
    mub <- fam$linkinv(etab)
    mu.eta.val <- fam$mu.eta(etab)
    mf$yy <- etab + (fit0$y - mub)/mu.eta.val - off
    wz <- w * mu.eta.val^2/fam$variance(mub)	
    ########################################
    mf$wts <- wz
    mf$invwt <- 1 / wz
    mcall$correlation <- cor_tree(value = vc.init, paras = list(V=V, W=mf$invwt))
    
    mcall$data <- mf
  }
  attributes(fit$logLik) <- NULL
  fit$call <- Call
  fit$family <- family
  fit$logLik <- as.numeric(NA)
  fit$etab <- etab
  fit$hyper <- hyper
  fit$wts <- wz
  fit$yy <- yy
  fit$hest <- hest
  if (norm.kernel) {
    fit$best <- as.vector(V.tree %*% t(Zs * ZZ) %*% (solve(V2)) %*% (residuals(fit, type='response'))) 	# Shrinkaged b
  } else {
    fit$best <- as.vector(V.tree %*% t(ZZ) %*% (solve(V2)) %*% (residuals(fit, type='response'))) 	# Shrinkaged b
  }
  
  fit$Zm <- Zm
  fit$V.tree <- V.tree
  fit$lambda1 <- lambda1
  fit$lambda2 <- lambda2
  fit$lambda3 <- lambda3
  fit$converge <- TRUE
  fit$norm.kernel <- norm.kernel
  #oldClass(fit) <- c("glmmTree", "glmmPQL", oldClass(fit))
  fit
}




glmmTreeb <-function(y,Z,X=NULL,D,lambda1,lambda2){
  
  nlambda1.0=nlambda1=length(lambda1)
  nlambda2.0=nlambda2=length(lambda2)
  lambda1.0=lambda1
  lambda2.0=lambda2
  
  n=length(y)
  
  mind=matrix(1,nlambda1,nlambda2)
  dimnames(mind)=list(lambda1=lambda1,lambda2=lambda2)
  outs=list()
  
  for(i in 1:nlambda1){
    for(j in 1:nlambda2){
      if(is.null(X)){
        out=try(glmmTree(y ~ 1, family='binomial', lambda1=lambda1[i], lambda2=lambda2[j], D=D, Z=Z))
      }else{
        out=try(glmmTree(y ~ X, family='binomial', lambda1=lambda1[i], lambda2=lambda2[j], D=D, Z=Z))
      }
      if(length(out)==1){
        mind[i,j]=0
      }else{
        outs[[paste(lambda1[i])]][[paste(lambda2[j])]]<-out
      }
    }
  }
  
  if(sum(apply(mind,1,sum))==0 | sum(apply(mind,2,sum))==0) return(NULL)
  lambda1=lambda1[apply(mind,1,sum)!=0]
  lambda2=lambda2[apply(mind,2,sum)!=0]
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  mind=as.matrix(mind[apply(mind,1,sum)!=0,apply(mind,2,sum)!=0,drop=F])
  
  
  if(nlambda1>0){
    Zms=V.trees=betas=us=etabs=list()
    for(i in 1:nlambda1){
      for(j in 1:nlambda2){
        if(mind[i,j]==1){
          betas[[paste(lambda1[i])]][[paste(lambda2[j])]]=as.matrix(coef(outs[[paste(lambda1[i])]][[paste(lambda2[j])]]))
          us[[paste(lambda1[i])]][[paste(lambda2[j])]]=as.matrix(outs[[paste(lambda1[i])]][[paste(lambda2[j])]]$best)
          Zms[[paste(lambda1[i])]][[paste(lambda2[j])]]=as.matrix(outs[[paste(lambda1[i])]][[paste(lambda2[j])]]$Zm)
          V.trees[[paste(lambda1[i])]][[paste(lambda2[j])]]=as.matrix(outs[[paste(lambda1[i])]][[paste(lambda2[j])]]$V.tree)
          etabs[[paste(lambda1[i])]][[paste(lambda2[j])]]=as.matrix(outs[[paste(lambda1[i])]][[paste(lambda2[j])]]$etab)
        }
      }
    }
  }else{
      warning("No solution for current parameter setting!")
      return(NULL)
  }
  
  list(betas=betas,us=us,Zms=Zms,V.trees=V.trees,etabs=etabs,lambda1=lambda1,lambda2=lambda2,nlambda1=nlambda1,nlambda2=nlambda2,mind=mind,lambda1.0=lambda1.0,lambda2.0=lambda2.0)

}
  


predict.glmmTreeb<-function(obj,X=NULL,Z){
    
    lambda1=obj$lambda1
    lambda2=obj$lambda2
    nlambda1=length(lambda1)
    nlambda2=length(lambda2)
    mind=obj$mind
    etas=list()
    
    if (is.vector(Z)) {
      Z=matrix(Z, nrow=1)
    }
    
    if(is.null(X)){
      X=as.matrix(rep(1,nrow(Z)))
    }else{
      X=cbind(1,X)
    }
    
    if(sum(mind)==0){
      etas[[1]]=as.matrix(rep(1,nrow(Z)))
      return(etas)
    }
    
    for(i in 1:nlambda1){
      ZZ=Z^lambda1[i]
      ZZ[Z==0]=0
      eta=NULL
      for(j in 1:nlambda2){
        if(mind[i,j]==1){
          newZZ=sweep(ZZ, 2, obj$Zms[[paste(lambda1[i])]][[paste(lambda2[j])]],'-') 
          Zs= 1 / sqrt(diag(newZZ %*% obj$V.trees[[paste(lambda1[i])]][[paste(lambda2[j])]] %*% t(newZZ)))
          e=X %*% obj$betas[[paste(lambda1[i])]][[paste(lambda2[j])]] + (Zs * newZZ) %*% obj$us[[paste(lambda1[i])]][[paste(lambda2[j])]]
          e=exp(e)/(1+exp(e))
          eta=cbind(eta,e)
        }
      }
      colnames(eta)=lambda2[mind[i,]!=0]
      etas[[paste(lambda1[i])]]=eta
      
    }
    
    etas
    
}



getmin2 <- function(lambda,cvm,cvsd){
  cvmin <- min(cvm,na.rm=TRUE)
  idmin <- cvm<=cvmin
  lambda.min <- max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,cvmin=cvmin)
}


setcrit <- function(obj,X,Z,y,crit,whichi){
  pred <- predict.glmmTreeb(obj, X,Z)
  prob_min <- 1e-5
  prob_max <- 1-prob_min	
  for(i in seq(obj$nlambda1)) crit[whichi,i,match(colnames(pred[[i]]),obj$lambda2)] <-pmin(pmax(pred[[i]],prob_min),prob_max)
  crit
}


auc <- function(Y,predmat){
  if(is.null(ncol(predmat)))
    predmat <- as.matrix(predmat)
  rprobmat <- apply(predmat,2,rank)
  n1 <- sum(Y);n0 <- length(Y)-n1	
  if(n1==0 || n0==0 ) {
    warning("No two classes within fold.")
    return(rep(0,length(Y)))
  }
  R1 <- apply(rprobmat[Y==1,,drop=FALSE],2,sum)
  umat <- R1-n1*(n1+1)/2
  umat <- umat/(n1*n0) #ncol(predmat) is length of lambda1 for current lambda2
  auc <- pmax(umat,1-umat)
  return(auc)
}



setcritcv <- function(crit,cv.ind,lambda1,lambda2,nfolds,y){
  
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  outmat <-  array(NA,c(nfolds,nlambda1,nlambda2))
  good  <- array(0,c(nfolds,nlambda1,nlambda2))
  dimnames(outmat)=dimnames(good)=list(nfolds=1:nfolds,lambda1=lambda1,lambda2=lambda2)
  
  for(ifold in seq(nfolds)){
    whichi <- cv.ind==ifold
    for(i in seq(nlambda1)){  # different #lambda2 for different fold
      id=colSums(apply(crit[whichi,i,,drop=FALSE],2,is.na))==0
      if(sum(id)!=0){
        outmat[ifold,i,id]=auc(y[whichi],crit[whichi,i,id])
      }
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



predict.cv.glmmTreeb<-function(obj.cv,X,Z){
  predict.glmmTreeb(obj.cv$obj.min,X,Z)[[1]]
}


cv.glmmTreeb<-function(y,Z,X=NULL,D,lambda1,lambda2, nfolds=5,trace=TRUE){
  
  obj=glmmTreeb(y,Z,X,D,lambda1,lambda2)
  if (is.null(obj)){
    warning("No solution for current parameter setting!")
    return(NULL)
  } 
  n=length(y)
  lambda1=obj$lambda1
  lambda2=obj$lambda2
  nlambda1=length(lambda1)
  nlambda2=length(lambda2)
  mind=obj$mind
  
  crit <- array(NA, c(n,nlambda1, nlambda2), dimnames=list(y=1:n, paste(lambda1),paste(lambda2)))
  
  if (min(table(y)) < nfolds) stop("nfolds is larger than the smaller of 0/1 in the data; decrease nfolds")
  if( (n/nfolds <10) && type.measure=="auc"){
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.glmgraph; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    type.measure="deviance"
  }
  if( (n/nfolds <3)){
    warning("Option grouped=FALSE enforced in cv.glmgraph, since < 3 observations per fold",call.=FALSE)
  }
  ind1 <- which(y==1)
  ind0 <- which(y==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
  cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
  cv.ind <- numeric(n)
  cv.ind[y==1] <- cv.ind1
  cv.ind[y==0] <- cv.ind0
  
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
    
    obj.i <- do.call(glmmTreeb,cv.args)
    
    if(is.null(obj.i)) stop("No solution within nfold, please reset your folds again!")
    X2 <- X[whichi, , drop=FALSE]
    Z2 <- Z[whichi, , drop=FALSE]
    y2 <- y[whichi]
    crit=setcrit(obj.i,X2,Z2,y2,crit,whichi)
  }
  
  
  critcv=setcritcv (crit,cv.ind,lambda1,lambda2,nfolds,y)
  
  
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
    lmin <- getmin2(as.numeric(colnames(cvm[[i]])),-cvm[[i]],cvsd[[i]])
    cvmat[i,] <- c(lambda1[i],lmin$lambda.min, lmin$cvmin)
  } 
  
  if(sum(cvmat)==0){
    rowids=unique(which(critcv$Nmat==max(critcv$Nmat), arr.ind=TRUE)[,1])
    for(i in 1:length(rowids)){
      ind <- critcv$Nmat[rowids[i],] >=max(critcv$Nmat)
      cvm[[i]] <- apply(critcv$crit[,rowids[i],ind,drop=FALSE],c(2,3),mean,na.rm=TRUE)
      cvsd[[i]] <- sqrt( apply(scale(critcv$crit[,rowids[i],ind],cvm[[i]],FALSE)^2,2,mean,na.rm=TRUE)/critcv$Nmat[rowids[i],ind])
      lmin <- getmin2(as.numeric(colnames(cvm[[i]])),-cvm[[i]],cvsd[[i]])
      cvmat[rowids[i],] <- c(lambda1[rowids[i]],lmin$lambda.min, lmin$cvmin)
    }
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
  obj.min <- do.call(glmmTreeb, cv.args)
  
  obj.cv <- list(obj.min=obj.min,cvmat=cvmat,cvmin=cvmin,cvm=cvm,cvsd=cvsd,lambda1.min=lambda1.min,lambda2.min=lambda2.min,
                 beta.min=obj.min$betas[[1]],u.min=obj.min$us[[1]])
  obj.cv
  
  
}







