######=====================================================================================
###### library of base functions
######=====================================================================================


gls2 <- function (model, data = sys.frame(sys.parent()), correlation = NULL, 
		weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail, 
		control = list(), verbose = FALSE, lower= -Inf, upper = Inf) {
	
	# Copied from nlmb:gls, allowing more controls
	Call <- match.call()
	controlvals <- glsControl()
	if (!missing(control)) {
		controlvals[names(control)] <- control
	}
	if (!inherits(model, "formula") || length(model) != 3L) {
		stop("\nmodel must be a formula of the form \"resp ~ pred\"")
	}
	method <- match.arg(method)
	REML <- method == "REML"
	if (!is.null(correlation)) {
		groups <- getGroupsFormula(correlation)
	}
	else groups <- NULL
	glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
	model <- terms(model, data = data)
	mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
					groups), data = data, na.action = na.action)
	if (!missing(subset)) {
		mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
	}
	mfArgs$drop.unused.levels <- TRUE
	dataMod <- do.call("model.frame", mfArgs)
	origOrder <- row.names(dataMod)
	if (!is.null(groups)) {
		groups <- eval(parse(text = paste("~1", deparse(groups[[2L]]), 
								sep = "|")))
		grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups, 
								asList = TRUE)))
		ord <- order(grps)
		grps <- grps[ord]
		dataMod <- dataMod[ord, , drop = FALSE]
		revOrder <- match(origOrder, row.names(dataMod))
	}
	else grps <- NULL
	X <- model.frame(model, dataMod)
	contr <- lapply(X, function(el) if (inherits(el, "factor")) 
					contrasts(el))
	contr <- contr[!unlist(lapply(contr, is.null))]
	X <- model.matrix(model, X)
	if (ncol(X) == 0L) 
		stop("no coefficients to fit")
	y <- eval(model[[2L]], dataMod)
	N <- nrow(X)
	p <- ncol(X)
	parAssign <- attr(X, "assign")
	fTerms <- terms(as.formula(model), data = data)
	namTerms <- attr(fTerms, "term.labels")
	if (attr(fTerms, "intercept") > 0) {
		namTerms <- c("(Intercept)", namTerms)
	}
	namTerms <- factor(parAssign, labels = namTerms)
	parAssign <- split(order(parAssign), namTerms)
	attr(glsSt, "conLin") <- list(Xy = array(c(X, y), c(N, ncol(X) + 
									1L), list(row.names(dataMod), c(colnames(X), deparse(model[[2]])))), 
			dims = list(N = N, p = p, REML = as.integer(REML)), logLik = 0)
	glsEstControl <- controlvals["singular.ok"]
	glsSt <- Initialize(glsSt, dataMod, glsEstControl)
	parMap <- attr(glsSt, "pmap")
	numIter <- numIter0 <- 0L
	repeat {
		oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
		if (length(coef(glsSt))) {
			optRes <- if (controlvals$opt == "nlminb") {
						nlminb(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
											glsPars), control = list(trace = controlvals$msVerbose, 
										iter.max = controlvals$msMaxIter), lower = lower, upper = upper)
					}
					else {
						optim(c(coef(glsSt)), function(glsPars) -logLik(glsSt, 
											glsPars), method = controlvals$optimMethod, 
								control = list(trace = controlvals$msVerbose, 
										maxit = controlvals$msMaxIter, reltol = if (numIter == 
														0L) controlvals$msTol else 100 * .Machine$double.eps),
								lower = lower, upper = upper)
					}
			coef(glsSt) <- optRes$par
		}
		else {
			optRes <- list(convergence = 0)
		}
		attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
		if (!needUpdate(glsSt)) {
			if (optRes$convergence) 
				stop(optRes$message)
			break
		}
		numIter <- numIter + 1L
		glsSt <- update(glsSt, dataMod)
		aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
		conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
		aConv <- c(beta = max(conv[1:p]))
		conv <- conv[-(1:p)]
		for (i in names(glsSt)) {
			if (any(parMap[, i])) {
				aConv <- c(aConv, max(conv[parMap[, i]]))
				names(aConv)[length(aConv)] <- i
			}
		}
		if (verbose) {
			cat("\nIteration:", numIter)
			cat("\nObjective:", format(optRes$value), "\n")
			print(glsSt)
			cat("\nConvergence:\n")
			print(aConv)
		}
		if (max(aConv) <= controlvals$tolerance) {
			break
		}
		if (numIter > controlvals$maxIter) {
			stop("maximum number of iterations reached without convergence")
		}
	}
	glsFit <- attr(glsSt, "glsFit")
	namBeta <- names(glsFit$beta)
	attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
			glsFit$varBeta * sqrt((N - REML * p)/(N - p))
	varBeta <- crossprod(varBeta)
	dimnames(varBeta) <- list(namBeta, namBeta)
	Fitted <- fitted(glsSt)
	if (!is.null(grps)) {
		grps <- grps[revOrder]
		Fitted <- Fitted[revOrder]
		Resid <- y[revOrder] - Fitted
		attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
	}
	else {
		Resid <- y - Fitted
		attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
	}
	names(Resid) <- names(Fitted) <- origOrder
	if (controlvals$apVar) {
		apVar <- glsApVar(glsSt, glsFit$sigma, .relStep = controlvals[[".relStep"]], 
				minAbsPar = controlvals[["minAbsParApVar"]], natural = controlvals[["natural"]])
	}
	else {
		apVar <- "Approximate variance-covariance matrix not available"
	}
	dims <- attr(glsSt, "conLin")[["dims"]]
	dims[["p"]] <- p
	attr(glsSt, "conLin") <- NULL
	attr(glsSt, "glsFit") <- NULL
	estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
			coefficients = glsFit[["beta"]], varBeta = varBeta, sigma = glsFit$sigma, 
			apVar = apVar, logLik = glsFit$logLik, numIter = if (needUpdate(glsSt)) numIter else numIter0, 
			groups = grps, call = Call, method = method, fitted = Fitted, 
			residuals = Resid, parAssign = parAssign, na.action = attr(dataMod, 
					"na.action"))
	if (inherits(data, "groupedData")) {
		attr(estOut, "units") <- attr(data, "units")
		attr(estOut, "labels") <- attr(data, "labels")
	}
	attr(estOut, "namBetaFull") <- colnames(X)
	class(estOut) <- "gls"
	estOut
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

	diag(V) <- exp(hyper) * diag(V) + W	
#	V <- nearPD(V)$mat
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


glmmTree <- function (model, family, data, weights, 
		lambda1, lambda2, lambda3=2, D, Z, vc.init=0, norm.kernel=TRUE, lower = -20, upper = 20,
		control, niter = 10, verbose = TRUE, ...) {
	# The code is modified on the glmmPQL from MASS package
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
	
#	if (!missing(correlation) && !is.null(attr(correlation, "formula"))) 
#		allvars <- c(allvars, all.vars(attr(correlation, "formula")))
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
	keep <- is.element(nm, c("model", "data", "subset", "lower", "upper", 
					"na.action", "control"))
	for (i in nm[!keep]) mcall[[i]] <- NULL
	model[[2L]] <- quote(yy)
	mcall[["model"]] <- model
#	mcall[[1L]] <- quote(nlme::gls)
    mcall[[1L]] <- quote(gls2)
	mcall$method <- "ML"
	
#	mcall$weights <- quote(varFixed(~invwt))
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
			message(gettextf("iteration %d", i), domain = NA)
		err <- try(fit <- eval(mcall))
		if (class(err) == 'try-error') {
			fit <- list(converge = FALSE)
			return(fit)
		}
	    
		#######################################
		# Modified based on PQL
	    # hest: Microbiome effect
		hyper <- coef(fit$modelStruct$corStruct)
#		cat(hyper, '\n')
#		cat(fit$sigma)
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
	oldClass(fit) <- c("glmmTree", "glmmPQL", oldClass(fit))
	fit
}


predict.glmmTree <- function (object, newX = NULL, newZ = NULL, type = c('link', 'response')) {
	type <- match.arg(type)
    
	if (object$converge == TRUE) {
		if (missing(newZ)) {
			pred <- object$etab
			pred <- switch(type, link = pred, response = object$family$linkinv(pred))
		} else {
			if (is.vector(newZ)) {
				newZ <- matrix(newZ, nrow=1)
				newX <- matrix(c(1, newX), nrow=1)
			} else {
				newX <- cbind(rep(1, nrow(newZ)), newX)
			}
			temp <- newZ ^ object$lambda1
			temp[newZ == 0] <- 0
			newZ <- sweep(temp, 2, object$Zm, '-') 
			if (object$norm.kernel) {
				Zs <- 1 / sqrt(diag(newZ %*% object$V.tree %*% t(newZ)))
				pred <- as.vector(newX %*% matrix(coef(object)) + (Zs * newZ) %*% object$best)
			} else {
				pred <- as.vector(newX %*% matrix(coef(object)) + newZ %*% object$best)
			}

			switch(type, response = {
						pred <- object$family$linkinv(pred)
					}, link = pred)
		}
		pred
	} else {
		pred <- NA
	}
    pred
}

