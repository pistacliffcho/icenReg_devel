###			SEMIPARAMETRIC UTILITIES

adjustIntervals <- function(B = c(0,1), yMat, eps = 10^-10){
  isCensored <- yMat[,2] - yMat[,1] > (2 * eps)
  if(B[1] == 0) yMat[isCensored,1] = yMat[isCensored,1] + eps
  if(B[2] == 0) yMat[isCensored,2] = yMat[isCensored,2] - eps
  return(yMat)
}

findMaximalIntersections <- function(lower, upper){
	allVals <- sort(unique(c(lower,upper)) )
	isLeft <- allVals %in% lower
	isRight <- allVals %in% upper
	miList <- .Call("findMI", allVals, isLeft, isRight, lower, upper)
	names(miList) <- c('l_inds', 'r_inds', 'mi_l', 'mi_r')	
	return(miList)
}


bs_sampleData <- function(rawDataEnv, weights){
	n <- length(rawDataEnv[['y']][,1])
	sampEnv <- new.env()
	sampInds <- sample(1:n, ceiling(sum(weights)), replace = TRUE, prob = weights)
	tabledInds <- table(sampInds)
	unqInds <- as.numeric(names(tabledInds))
	weights <- as.numeric(tabledInds)
	sampEnv[['x']] <- rawDataEnv[['x']][unqInds,]
	sampEnv[['y']] <- rawDataEnv[['y']][unqInds,]
	sampEnv[['w']] <- weights
	return(sampEnv)
}

getBS_coef <- function(sampDataEnv, callText = 'ic_ph', other_info){ 
	xMat <- cbind(sampDataEnv$x,1)
	invertResult <- try(diag(solve(t(xMat) %*% xMat )), silent = TRUE)
	if(is(invertResult, 'try-error')) {return( rep(NA, ncol(xMat) -1) ) }
	output <- fit_ICPH(sampDataEnv$y, sampDataEnv$x, callText, sampDataEnv$w, other_info)$coefficients
	return(output)
}


expandX <- function(formula, data, fit){
  if(inherits(fit, 'ic_np')) return(NULL)
	tt <- terms(fit)
	Terms <- delete.response(tt)
	 m <- model.frame(Terms, as.data.frame(data), na.action = na.pass, xlev = fit$xlevels)
	 x <- model.matrix(Terms, m)
	 if(length(x) == 0) return(NULL)
	 ans <- as.matrix(x[,-1])
	 if(nrow(ans) != nrow(x)){
	   ans <- t(ans)
	   if(nrow(ans) != nrow(x) ){stop('nrow(ans) != nrow(x)')}
	   if(ncol(ans) != (ncol(x) - 1) ) {stop('ncol(ans) != ncol(x) - 1')}
	 }
	 return(ans)
}

removeSurvFromFormula <- function(formula, ind = 2){ 
  if(formula[[ind]][[1]] == 'Surv'){
    if(formula[[ind]][[4]] != 'interval2') stop('Surv type not "interval2" in formula')
    formula[[ind]][[1]] <- as.name('cbind')
    formula[[ind]][[4]] <- NULL
  }
  return(formula)
}

expandY <- function(formula, data, fit){
  if(inherits(fit, 'ic_np')){
    if(ncol(data) != 2) stop('expandY expected an nx2 matrix for data')
    return(data)
  }
  if(is.null(data)) data <- fit$getRawData()
  newFormula <- formula
  newFormula[[3]] <- newFormula[[2]]
  newFormula[[2]] <- NULL
  newFormula <- removeSurvFromFormula(newFormula, 2)
  ans <- model.matrix(newFormula, data)
  ans <- ans[,-1] #Dropping Intercept
  return(ans)
}

getResponse <- function(fit, newdata = NULL){
  if(is.null(newdata))
    newdata = fit$getRawData()
  ans <- expandY(fit$formula, newdata, fit)
  return(ans)
}
###		PARAMETRIC FIT UTILITIES

make_par_fitList <- function(y_mat, x_mat, parFam = "gamma", 
                             link = "po", leftCen = 0, rightCen = Inf,
                             uncenTol = 10^-6, regnames,
                             weights, callText){
  k_reg <- getNumCovars(x_mat)
  etaOffset = 0
  if(!is.matrix(x_mat))
    x_mat <- matrix(x_mat, ncol = 1)
#  if(recenterCovar == TRUE){
#    prcomp_xmat <- prcomp(x_mat, center = TRUE, scale. = TRUE)
#    x_mat <- prcomp_xmat$x
#  }
  
  isUncen <- abs(y_mat[,2] - y_mat[,1]) < uncenTol
  mean_uncen_t <- (y_mat[isUncen,1] + y_mat[isUncen,2])/2
  y_mat[isUncen,1] <- mean_uncen_t
  y_mat[isUncen,2] <- mean_uncen_t
  isRightCen <- y_mat[,2] == rightCen
  isLeftCen <- y_mat[,1]  <= leftCen
  isGCen <- !(isUncen | isRightCen | isLeftCen)
  
  if(any(isRightCen & isUncen))	stop('estimator not defined if left side of interval = 0')
  if(any(isLeftCen & isUncen))	stop('uncensored times cannot be equal = 0. Try replacing exact times = 0 with really small numbers')
  if(any(y_mat[,1] > y_mat[,2]))	stop('left side of interval cannot be larger than right!')
  
  s_t <- unique( as.numeric(y_mat) )
  uncenInd_s <- match(y_mat[isUncen, 1], s_t)
  d_t <- unique(s_t[uncenInd_s])
  uncenInd_d <- match(y_mat[isUncen, 1], d_t)
  uncenInd_mat <- as.matrix(cbind(uncenInd_d, uncenInd_s))
  
  rightCenInd <- match(y_mat[isRightCen,1], s_t)
  leftCenInd <- match(y_mat[isLeftCen,2], s_t)
  
  leftSideInd <- match(y_mat[isGCen, 1], s_t)
  rightSideInd <- match(y_mat[isGCen, 2], s_t)
  
  gicInd_mat <- as.matrix(cbind(leftSideInd, rightSideInd))
  
  w_reordered <- c(weights[isUncen], weights[isGCen], weights[isLeftCen], weights[isRightCen])
  
  if(is.matrix(x_mat)	){
    if(ncol(x_mat) > 1)	x_mat_rearranged <- rbind(x_mat[isUncen,], x_mat[isGCen,], x_mat[isLeftCen,], x_mat[isRightCen,])
    else				x_mat_rearranged <- matrix(c(x_mat[isUncen], x_mat[isGCen], x_mat[isLeftCen], x_mat[isRightCen]), ncol = 1)
  }
  else if(length(x_mat) != 0)		x_mat_rearranged <- matrix(c(x_mat[isUncen], x_mat[isGCen], x_mat[isLeftCen], x_mat[isRightCen]), ncol = 1)
  else							x_mat_rearranged <- matrix(ncol = 0, nrow = nrow(x_mat))
  storage.mode(x_mat_rearranged) <- 'double'
  x_mat_rearranged <- as.matrix(x_mat_rearranged)	
  
  if(k_reg == 0)	x_mat_rearranged <- matrix(nrow = nrow(x_mat), ncol = k_reg)
  
  #regnames = colnames(x_mat_rearranged)
  if(parFam == 'gamma') {parInd = as.integer(1); k_base = 2; bnames = c('log_shape', 'log_scale')}
  else if(parFam == 'weibull') {parInd = as.integer(2); k_base = 2; bnames = c('log_shape', 'log_scale')}
  else if(parFam == 'lnorm') {parInd = as.integer(3); k_base = 2; bnames = c('mu', 'log_s')}
  else if(parFam == 'exponential') {parInd = as.integer(4); k_base = 1; bnames = 'log_scale'}
  else if(parFam == 'loglogistic') {parInd = as.integer(5); k_base = 2; bnames = c('log_alpha', 'log_beta')}
  else if(parFam == 'generalgamma') {parInd = as.integer(6); k_base = 3; bnames = c('mu', 'log_s', 'Q')}
  else stop('parametric family not supported')
  
  hessnames = c(bnames, regnames)
  
  if(link == 'po') linkInd = as.integer(1)
  else if (link == 'ph') linkInd = as.integer(2)
  else if (link == 'aft') linkInd = as.integer(3)
  else stop('link function not supported')
  
  hessian <- matrix(numeric(), nrow = (k_reg + k_base), ncol = (k_reg + k_base))

  
  ans <- list(
    s_t               = s_t,
    d_t               = d_t,
    covars            = x_mat_rearranged,
    uncenInd_mat      = uncenInd_mat,
    gicInd_mat        = gicInd_mat,
    leftCenInd        = leftCenInd,
    rightCenInd       = rightCenInd,
    parInd            = parInd,
    linkType          = linkInd,
    hessian           = hessian,
    w                 = as.numeric(w_reordered),
    bnames            = bnames,
    regnames          = regnames,
    hessnames         = hessnames
  )
  
  return(ans)
}



###			IMPUTATION UTILITIES

imputeCensoredData_exp <- function(l, u, impInfo, dist = 'web', maxVal){
	
	if(dist == 'exp'){
		rate <- impInfo$indInters
		p_l <- pexp(l, rate)
		p_u <- pexp(u, rate)
		
		samp_q <- runif(length(l), p_l, p_u)
		samp_val <- qexp(samp_q, rate)
		is.inf <- samp_val == Inf
		samp_val[is.inf] <- u[is.inf]
		return(samp_val)
	}
	if(dist =='web'){
		inter <- impInfo$indInters
		scale <- impInfo$scale
		p_l <- pweibull(l, shape = 1/scale, scale = inter)
		p_u <- pweibull(u, shape = 1/scale, scale = inter)
		samp_q <- runif(length(l), p_l, p_u)
		samp_val <- qweibull(samp_q, shape = 1/scale, scale = inter)
		is.inf <- samp_val == Inf
		samp_val[is.inf] <- u[is.inf]
		replaceInds <- samp_val > maxVal
		samp_val[replaceInds] <- maxVal
		return(samp_val)
	}
}


fullParamFit <- function(formula, data, param_y, dist = 'weibull'){
	data$param_y = param_y
	fit <- survreg(formula, data, dist = dist, x = TRUE)
	fit$var_chol <- chol(fit$var)
	return(fit)
}

fullParamFit_exp <- function(formula, data, param_y, rightCenVal, dist = 'weibull'){
	data$param_y = param_y
	isInf <- data$param_y[,2] == Inf
	data$param_y[isInf,2] <- rightCenVal
	if( any(data$param_y[,2] > rightCenVal) )	stop('finite right side of interval greater than provided rightCenVal. Use larger rightCenVal')
	
	fit <- survreg(formula, data, dist = dist, x = TRUE)
	fit$var_chol <- chol(fit$var)
	return(fit)
}

simPars_fromFit <- function(fit, web = TRUE){
	means <- fit$coef
	if(web)
		means = c(means, log(fit$scale) )
	
	k <- length(means)
	sampPars <- means + fit$var_chol %*% rnorm(k)
	if(web){
		scale <- exp(sampPars[k])
		sampPars <- sampPars[-k]
	}
	indInters <- exp(fit$x %*% sampPars)
	output <- list(indInters = indInters)
	if(web)
		output$scale = scale
	return(output)
}





#####		SIMULATION UTILITIES


simRawTimes <- function(b1 = 0.5, b2 = -0.5, n = 100, shape1 = 2, shape2 = 2){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ^(1/nu)
    trueTimes <- qbeta(adjQ, shape1 = shape1, shape2 = shape2)
	return(data.frame(y = trueTimes, x1 = x1, x2 = x2, obs = rep(1, n)))
}

simRawExpTimes <- function(b1 = 0.5, b2 = -0.5, n = 100, rate = 1){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ^(1/nu)
    trueTimes <- qexp(adjQ, rate = rate)
	return(data.frame(y = trueTimes, x1 = x1, x2 = x2, obs = rep(1, n)))
}




####		DIAGNOSTIC UTILITIES
subSampleData <- function(data, max_n_use, weights){
	if(nrow(data) <= max_n_use) return(list( data = data, w = weights))
	if(is.null(weights)) weights <- rep(1, nrow(data))
	sampInd <- sample(1:nrow(data), max_n_use, prob = weights)
	tabledInds <- table(sampInd)
	
	newWeights <- as.numeric(names(tabledInds))
	newInds <- as.numeric(tabledInds)
	subData <- data[newInds,]
	return(list(data = subData, w = newWeights))
}

removeVarFromCall <- function(call, varName){
	newcall <- call
	charCall <- as.character(newcall)
	varName_plus <- paste(varName, '+')
	plus_varName <- paste('+', varName)
	varName_times <- paste('*', varName)
	times_varName <- paste(varName, '*')
	charCall[3] <- gsub(varName_plus, '', charCall[3], fixed = TRUE)
	charCall[3] <- gsub(plus_varName, '', charCall[3], fixed = TRUE)
	charCall[3] <- gsub(varName_times, '', charCall[3], fixed = TRUE)
	charCall[3] <- gsub(times_varName, '', charCall[3], fixed = TRUE)
	charCall[3] <- gsub(varName, '', charCall[3])
	if(nchar(charCall[3]) == 0)	charCall[3] <- '0'
	pastedCall <- paste(charCall[2], charCall[1], charCall[3])
	newcall <- try(parse(text = pastedCall)[[1]], silent = TRUE)
	if(inherits(newcall, 'try-error')){
		warnMessage <- paste('Unable to determine formula for splitting on variable  ', varName, 
		#		'Original formula = ', call,
				'\nAttempted formula after removing variable = ', pastedCall,
				'\nPossibly due to interaction term. Currently splitting on an interaction term is not supported. User will have to split and diagnosis manually')
				stop(warnMessage)
	}
	newcall
}

splitData <- function(data,		#full data
					 varName,	#name to split data on
					 splits,	#list of information necessary to make splits
					 splitFun,	#function that uses splits elements to divide up data
					 weights
					 ){
					 splitData <- new.env()
					 splitNames <- names(splits)
					 for(sn in splitNames)
					 	splitData[[sn]] <- splitFun(sn, varName, data, splits, weights)	
	return(splitData)
}

makeFactorSplitInfo <- function(vals, levels){
	sInfo <- list()
	sInfo$splits <- list()
	sInfo$splits[levels] <- 0
	sInfo$splitFun <- function(splitName, varName, data, splits, weights){
		keepInds <- which(data[[varName]] == splitName) 
		return(list(data = data[keepInds,], w = weights[keepInds]) )
	}
	sInfo
}

makeNumericSplitInfo <- function(vals, cuts){
	sInfo <- list()
	sInfo$splits <- list()
	prep_cuts <- unique(sort(c(-Inf, cuts, Inf)))
	for(i in 1:(length(prep_cuts) - 1) ){
		theseCuts <- prep_cuts[c(i, i+1)]
		cutNames <- paste0('(', round(theseCuts[1], 2), ',', round(theseCuts[2],2), ']')
		sInfo$splits[[cutNames]] <- theseCuts
	}
	sInfo$splitFun <- function(splitName, varName, data, splits, weights){
		theseCuts <- splits[[splitName]]
		keep <- data[[varName]] > theseCuts[1] & data[[varName]] <= theseCuts[2]
		return(list(data = data[keep,], w = weights[keep]) )
	}
	sInfo
}

splitAndFit <- function(newcall, data, varName, splitInfo, fitFunction, model, weights){
	split_data <- splitData(data, varName = varName, splits = splitInfo$splits, 
							splitFun = splitInfo$splitFun, weights = weights)
	splitNames <- ls(split_data)
	splitFits <- new.env()
	for(sn in splitNames){
		theseData <- split_data[[sn]]$data
		theseWeights <- split_data[[sn]]$w
		splitFits[[sn]] <- ic_sp(newcall, data = theseData, model = model, weights = theseWeights)
	}
	return(splitFits)
}

s_exp <- function(x, par){return(1 - pexp(x, exp(-par[1]))) }
s_weib <- function(x, par){return(1 - pweibull(x, exp(par[1]), exp(par[2]))) }
s_gamma <- function(x, par){return(1 - pgamma(x, shape = exp(par[1]), scale = exp(par[2]) ) )}
s_lnorm <- function(x, par){return(1 - pnorm(log(x), mean = par[1], sd = exp(par[2])))}
s_loglgst <- function(x, par){
	a <- exp(par[1])
	b <- exp(par[2])
	ans <- 1 - 1 /(1 + (x / a)^(-b) )
	return(ans)
}

get_etas <- function(fit, newdata = NULL, reg_pars = NULL){
  if(fit$par == 'non-parametric'){ans <- 1; names(ans) <- 'baseline'; return(ans)}
  if(is.null(reg_pars)){ reg_pars <- default_reg_pars(fit) }
  if(!is.matrix(reg_pars)) reg_pars <- matrix(reg_pars, ncol = length(default_reg_pars(fit)))
  if(is.null(newdata)){
    ans <- 1
	  names(ans) <- 'baseline'
	  return(ans)
	}
	if(identical(newdata, 'midValues')){
  	ans <- 1
  	names(ans) <- 'Mean Covariate Values'
  	return(ans)
	}
  if(identical(rownames(newdata), NULL) ) {rownames(newdata) <- as.character(1:icr_nrow(newdata))}
	grpNames <- rownames(newdata)
	reducFormula <- fit$formula
	reducFormula[[2]] <- NULL
	new_x <- expandX(reducFormula, newdata, fit)
	if(is.null(new_x)){
	  if(is.matrix(reg_pars))
  	  ans <- rep(1, nrow(reg_pars))
	  else
	    ans <- rep(1, length(reg_pars))
	  return(ans)
	}
	covarOffset <- fit$covarOffset
#	new_x <- new_x - covarOffset
	new_x <- subtractOffset(new_x, covarOffset)
	log_etas <- as.numeric( new_x %*% t(reg_pars) ) 	
	etas <- exp(log_etas)
	if(length(grpNames) == length(etas)){ names(etas) <- grpNames }
	return(etas)
}

get_s_fun <- function(fit){
	if(inherits(fit, 'par_fit')){
		par_fam = fit$par
		if(par_fam == 'exponential') return(s_exp)
		if(par_fam == 'weibull') return(s_weib)
		if(par_fam == 'gamma') return(s_gamma)
		if(par_fam == 'lnorm') return(s_lnorm)
		if(par_fam == 'loglogistic') return(s_loglgst)
	}
}

po_link <- function(s, nu){ nu * s / (s * (nu - 1) + 1)	}
ph_link <- function(s, nu){ s^nu }
no_link <- function(s , nu){ s }
get_link_fun <- function(fit){
	if(fit$model == 'po') return(po_link)
	if(fit$model == 'ph')	return(ph_link)
  if(fit$model == 'none') return(no_link)
	stop('model type not recognized. Should be "ph", "po" or "none"')
}


findUpperBound <- function(val = 1, x, s_fun, link_fun, fit, eta, baseline = NULL){
  if(is.null(baseline)){ baseline <- default_baseline(fit) }
	fval <- 1 - link_fun(s_fun(val, baseline), eta)
	tries = 0
	while(tries < 100 & fval < x){
		tries = tries + 1
		val <- val * 10
		fval <- 1 - link_fun(s_fun(val, baseline), eta)
	}
	if(fval < x)	stop('finding upper bound for quantile failed!')
	return(val)
}


subsetData_ifNeeded <- function(i, data){
	if(is.null(data)) return(NULL)
	if(ncol(data) > 1) return(data[i,])
	newdata <- data.frame(data[i,])
	colnames(newdata) <- colnames(data)
	return(newdata)
}


getVarNames_fromFormula <- function(formula)
	getVar_fromRHS(formula[[3]], NULL)

getVar_fromRHS <- function(rhs, names){
	if(length(rhs) == 1 ){
		if(rhs != '+' & rhs != '*' & rhs != ":")
			return( c(names, as.character(rhs)))
		return(names)
	}
	names <- NULL
	for(i in seq_along(rhs))
		names <- getVar_fromRHS(rhs[[i]], names)
	return(names)
}


getFormula <- function(object){
	if(inherits(object, 'formula'))					return(object)
	else if(inherits(object, 'icenReg_fit'))		return(object$formula)
	stop('object should be either a formula or an icenReg_fit object')
}

getData <- function(fit){
	ans <- fit$getRawData()
	if(is.null(ans)) stop('Could not find data from fit. Original model must be built with data argument (rather than variables found in the Global Environment) supplied to be retreivable')
	return(ans)
}

get_tbull_mid_q <- function(p, s_t, tbulls){
  if(any(is.na(c(p, s_t, tbulls)))) stop('NAs provided to get_tbull_mid_q')
	x_low <- tbulls[,1]
	x_hi <- tbulls[,2]
	x_range <- range(tbulls)
	k <- length(tbulls)
	
	ans <- numeric()
	
	for(i in seq_along(p)){
		this_p <- p[i]
		if(this_p < 0 | this_p > 1)	stop('attempting to find survival quantile for invalid probability')	
		if(this_p == 0)			ans[i] <- x_range[1]
		else if(this_p == 1) 	ans[i] <- x_range[2]
		else{
			this_s <- 1 - this_p
			s_ind <- min(which(s_t <= this_s))
	
			s_res <- this_s - s_t[s_ind]
			s_jump <- s_t[s_ind- 1] - s_t[s_ind]
			
			this_q_low <- x_low[s_ind]
			this_q_hi  <- x_hi[s_ind]
			if(s_jump > 0)
				ans[i] <- this_q_hi + (this_q_low - this_q_hi) * s_res/s_jump
			else
				ans[i] <- this_q_hi		
		}
	}
	return(ans)
}

get_tbull_mid_p <- function(q, s_t, tbulls){
	x_low <- tbulls[,1]
	x_hi <- tbulls[,2]
	x_range <- range(tbulls)
	k <- length(tbulls)
	
	ans <- numeric()
	
	max_x_low 	<- max(x_low)
	min_x_hi 	<- min(x_hi)
		
	for(i in 1:length(q)){
		this_q <- q[i]
		if(this_q <= x_range[1]) 			ans[i] <- 0
		else if(this_q > x_range[2])		ans[i] <- 1
		else{
			l_ind <- max(which(x_low <= this_q))
			if(this_q < min_x_hi)	u_ind <- 1
			else					u_ind <- max(which(x_hi <= this_q))
			
			if(l_ind == u_ind){
						ans[i] = 1 - s_t[l_ind]
						}
			else{
				p_jump <- s_t[u_ind] - s_t[l_ind]
				int_lng <- x_hi[l_ind] - x_low[u_ind]
				q_res <- this_q - x_low[u_ind]
								
				if(int_lng > 0)
					ans[i] <- 1 - (s_t[u_ind] - p_jump * q_res/int_lng)	
				else
					ans[i] <- 1 - (s_t[u_ind] - p_jump)
												
				if(ans[i] < 0 | ans[i] > 1) cat('error occurred ')								
														
			}
		}
	}
	return(ans)
}

PCAFit2OrgParFit <- function(PCA_info, PCA_Hessian, PCA_parEsts, numIdPars){
	tot_k = numIdPars + length(PCA_info$scale)
	if(tot_k != length(PCA_parEsts) )		stop('incorrect dimensions for PCAFit2OrgParFit')
	if(numIdPars == 0)		pcaTransMat <- PCA_info$rotation
	else{
		pcaTransMat <- diag(1, tot_k)
		pcaTransMat[(-1:-numIdPars), (-1:-numIdPars)] <- PCA_info$rotation
	}
	
	for(i in seq_along(PCA_info$scale))
		pcaTransMat[i + numIdPars, ] <- pcaTransMat[i + numIdPars,]/PCA_info$scale[i]
	ans <- list()
	ans$pars <- pcaTransMat %*% PCA_parEsts
	ans$var  <- pcaTransMat %*% -solve(PCA_Hessian) %*% t(pcaTransMat)
	return(ans)
}




getNumCovars <- function(object){
  dimAns <- dim(object)
  if(is.null(dimAns)){
    if(length(object) > 0 ) return(1)
    return(0)
  }
  if(length(dimAns) == 2) return(dimAns[2])
  stop('problem with getNumCovars')
}



regmod2int <- new.env()
regmod2int[['ph']] <- as.integer(1)
regmod2int[['po']] <- as.integer(2)
regmod2int[['aft']] <- as.integer(2) 
# Above this is not a typo!
# This is because aft models are handled fundamentially differently than 
# po and ph models. So it actually is being given all 1's for etas to keep the machinery
# built for po and ph in place, but not actually doing anything
regmod2int[['none']] <- as.integer(0)

basemod2int <- new.env()
basemod2int[['sp']] <- as.integer(0)
basemod2int[['gamma']] <- as.integer(1)
basemod2int[['weibull']] <- as.integer(2)
basemod2int[['weib']] <- as.integer(2)
basemod2int[['lnorm']] <- as.integer(3)
basemod2int[['exponential']] <- as.integer(4)
basemod2int[['loglogistic']] <- as.integer(5)
basemod2int[['generalgamma']] <- as.integer(6)

getSurvProbs <- function(times, etas, baselineInfo, regMod, baseMod){
  regInt <- regmod2int[[regMod]]
  if(is.null(regInt)) stop('regMod type not recognized')
  baseInt <- basemod2int[[baseMod]]
  if(is.null(baseInt)) stop('baseMod type not recognized')
  ans <- .Call('s_regTrans', as.double(times), as.double(etas), baselineInfo, regInt, baseInt)
  return(ans)
}

getSurvTimes <- function(p, etas, baselineInfo, regMod, baseMod){
  if(any(p < 0))
    stop('probabilities provided to getSurvTimes are less than 0')
  if(any(p > 1))
    stop('probabilities provided to getSurvTimes are greater than 1')
  regInt <- regmod2int[[regMod]]
  if(is.null(regInt)) stop('regMod type not recognized')
  baseInt <- basemod2int[[baseMod]]
  if(is.null(baseInt)) stop('baseMod type not recognized')
  ans <- .Call('q_regTrans', as.double(p), as.double(etas), baselineInfo, regInt, baseInt)
  return(ans)
}


getSamplablePars <- function(fit){
  if(is(fit, 'sp_fit') | is(fit, 'par_fit') | is(fit, 'bayes_fit')) 
    return(fit$coefficients)
}

getSamplableVar <- function(fit){
  ans <- fit$var
  if(is.null(ans))  stop('coefficient variance not found for fit. If ic_ph model was fit, make sure to set bs_samples > 100 to get bootstrap sample of variance')
  return(ans)
}

sampBayesPar <- function(fit){
  nRow <- nrow(fit$samples)
  samp_ind <- round(runif(n = 1, min = 0.5, max = nRow + 0.5) ) 
  ans <- fit$samples[samp_ind,]
  return(ans)
}

sampPars <- function(mean, var){
  chol_var <- chol(var)
  k <- length(mean)
  ans <- rnorm(k) %*% chol_var + mean
  return(ans)
}

getBSParSample <- function(fit){
  if(inherits(fit, 'ic_np')) return(NULL)
  nBS_samps <- nrow(fit$bsMat)
  if(is.null(nBS_samps) ){
    stop('no bootstrap samples generated so cannot sample parameter values')
  }
  thisInd <- sample(1:nBS_samps, 1)
  return(fit$bsMat[thisInd, ] )
}

setSamplablePars <- function(fit, coefs){
  if(!inherits(fit, 'ic_np')){ 
    fit$coefficients <- coefs
    if(!inherits(fit, 'sp_fit')){
      n_base <- length(fit$baseline)
      fit$baseline[1:n_base] <- coefs[1:n_base]
    }
    else n_base = 0
    if(length(coefs) > n_base) fit$reg_pars[1:(length(coefs) - n_base)] <- coefs[(n_base+1):length(coefs)]
  }
}


fastNumericInsert <- function(newVals, target, indices){
  if(storage.mode(newVals) != 'double') storage.mode(newVals) <- 'double'
  if(storage.mode(target) != 'double') stop('target of fastNumericInsert MUST have storage.mode = "double"')
  if(storage.mode(indices) != 'integer') storage.mode(indices) <- 'integer'
  
  invisible(.Call('fastNumericInsert', newVals, target, indices) )
}

fastMatrixInsert <- function(newVals, targMat, rowNum = NULL, colNum = NULL){
  if(is.null(colNum)){
    if(is.null(rowNum)) stop('need either rowNum or colNum')
    newIndices <- (1:length(newVals)-1) * nrow(targMat) + rowNum 
    return(fastNumericInsert(newVals, targMat, newIndices))
  }
  newIndices <- 1:length(newVals) + (colNum - 1) * nrow(targMat)
  fastNumericInsert(newVals, targMat, newIndices)
}

updateDistPars <- function(vals, max_n){
  vals <- as.numeric(vals)
  n_v <- length(vals)
  if(n_v == 1){
    return(rep(vals, max_n))
  }
  if(n_v != max_n) stop('parameters are not of equal length (or 1)')
  return(vals)
}

getMaxLength <- function(thisList){
  n <- 0
  for(i in seq_along(thisList)) n <- max(c(n, length(thisList[[i]] ) ) )
  return(n)
}


### PREPROCESSING TOOLS
addIfMissing <- function(val, name, list){
  if(is.null(list[[name]])) list[[name]] <- val
  return(list)
}

addListIfMissing <- function(listFrom, listInto){
  listFromNames <- names(listFrom)
  for(n in listFromNames) listInto <- addIfMissing(listFrom[[n]], n , listInto)
  return(listInto)
}

readingCall <- function(mf){
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$formula = quote(formula)
  mf$data = quote(data)
  mf$na.action = quote(na.pass)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  ans <- list(mt = mt, mf = mf)
  return(ans)
}

makeIntervals <- function(y, mf){
  yMat <- as.matrix(y)[,1:2]
  if(is(y, 'Surv')){
    rightCens <- mf[,1][,3] == 0
    yMat[rightCens,2] <- Inf
    exact <- mf[,1][,3] == 1
    yMat[exact, 2] = yMat[exact, 1]
  }
  else{
    rightIsNA <- is.na(yMat[,2])
    yMat[rightIsNA,2] <- Inf
  }
  storage.mode(yMat) <- 'double'
  return(yMat)
}

checkWeights <- function(weights, yMat){
  if(is.null(weights)) 				weights = rep(1, nrow(yMat))
  if(length(weights) != nrow(yMat))	stop('weights improper length')
  if(any(is.na(weights) > 0) )		stop('NAs not allowed in weights')
  if(any(weights < 0)	)				stop('negative weights not allowed')
  return(weights)
}

checkMatrix <- function(x){
  testMat <- cbind(x, 1)
  invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
  if(is(invertResult, 'try-error'))
    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
}

makeQQFit <- function(fit){
  if(!is(fit, 'par_fit')){stop("icenReg's QQ plots only for models fit by ic_par")}
  thisData <- fit$getRawData()
  if(fit$model == 'aft'){
    theseEtas <- get_etas(fit, thisData)
    baseData <- getResponse(fit)
    baseData[,1] <- baseData[,1] / theseEtas
    baseData[,2] <- baseData[,2] / theseEtas
    spfit <- ic_np(baseData)
  }
  else{
    thisFormula <- fit$formula
    thisModel <- fit$model
    theseRegPars <- fit$reg_pars
    theseControls <- makeCtrls_icsp(regStart = theseRegPars)
    theseControls$updateReg = F
    spfit <- ic_sp(formula = thisFormula, model = thisModel, data = thisData, controls = theseControls)
  }
  return(spfit)
}




formula_has_plus <- function(form){
  form_char <- as.character(form)
  return("+" %in% form_char)
}

getFactorFromData <- function(formExpression, data){
  form_char <- as.character(formExpression)
  if(form_char[1] == 'factor') form_char = form_char[2]
  ans <- as.factor(data[[form_char]])
  return(ans)
}

icr_nrow <- function(x){
  ans <- nrow(x)
  if(is.null(ans)){
    if(length(x) > 0) ans <- 1
  }
  return(ans)
}

icr_ncol <- function(x){
  ans <- ncol(x)
  if(is.null(ans)){
    if(length(x) > 0) ans <- 1
  }
  return(ans)
}

icr_colMeans <- function(x){
  if(icr_ncol(x) > 1) return( colMeans(x) )
  return(mean(x))
}

get_dataframe_row <- function(df, row){
  if(is.null(df)) return(NULL)
  if(ncol(df) == 0) return(NULL)
  if(ncol(df)  > 1) return(df[row, ])
  ans <- as.data.frame( df[row,] )
  colnames(ans) <- colnames(df)
  rownames(ans) <- rownames(df)[row]
  return(ans)
}


default_baseline <- function(fit){
  if(!inherits(fit, 'fit_bayes')) return(fit$baseline)
  return(fit$MAP_baseline)
}
default_reg_pars <- function(fit){
  if(!inherits(fit, 'fit_bayes')) return(fit$reg_pars)
  return(fit$MAP_reg_pars)
}


sample_in_interval <- function(fit, newdata, lower_time, upper_time){
  p_l <- getFitEsts(fit, newdata, q = lower_time)
  p_u <- getFitEsts(fit, newdata, q = upper_time)
  raw_p <- runif(length(p_l), min = p_l, max = p_u)
  ans <- getFitEsts(fit, newdata, p = raw_p)
  return(ans)
}

sample_pars <- function(fit, samples = 100){
  if(is(fit, 'par_fit') | is(fit, 'sp_fit')){
    chol_var <- chol(vcov.icenReg_fit(fit))
    mean_coefs <- as.numeric( fit$coefficients )
    k <- length(mean_coefs)
    norm_samps <- matrix(rnorm(samples * k), nrow = samples)
    ans <- norm_samps %*% chol_var + rep(1, samples) %*% t(mean_coefs)
    return(ans)
  }
  if(is(fit, 'bayes_fit')){
    n_samps <- nrow(fit$samples)
    resamp = samples > n_samps
    samp_inds <- sample(1:n_samps, samples, replace = resamp)
    ans <- fit$samples[samp_inds,]
    return(ans)
  }
}

sample_etas_and_base <- function(fit, samples, newdata){
  all_pars <- sample_pars(fit, samples)
  nBase <- length(fit$baseline)
  basePars <- all_pars[,1:nBase]
  if(nBase == 1) basePars <- as.matrix(basePars, ncol = 1)
  regPars  <- all_pars[,-(1:nBase)]
  etas <- get_etas(fit, newdata, reg_pars = regPars)
  if(length(etas) == 0) etas <- rep(1, samples)
  ans <- list(baseMat = basePars, 
              etas = etas)
  return(ans)
}

subtractOffset <- function(new_x, offset){
  ncol_x <- ncol(new_x)
  ncol_off <- ncol(offset)
  if(ncol_x != ncol_off) stop('ncol(new_x) != ncol(offset)')
  for(i in seq_along(offset)){
    new_x[,i] <- new_x[,i] - offset[1,i]
  }
  return(new_x)
}


# This function throws an error if a user tries to include "cluster(x)" 
# on the right hand side of a formula
checkFor_cluster = function(form){
  # Extracting right hand side
  rhs = form[[3]]
  # Turning into characters
  rhs_char = as.character(rhs)
  
  if(any(grepl("cluster\\(", rhs_char) ) ) 
    stop("cluster(covar) not implemented in icenReg. To account for repeated measures, see ?ir_clustBoot")
}


# This function prepares response + feature matrices
#' @param frml formula object
#' @param data data.frame
make_xy = function(frml, df){
  ans = list()
  # Making a model.frame
  mod_frame = model.frame(formula = eval(frml), 
                          data = as.data.frame(df), 
                          # NA's allowed in y (right censoring) but not x
                          # We check x for nas manually
                          na.action = na.pass)
  # Getting x from model frame
  # trapping weird problem happens with diag_covar if rightside is ~0 
  x <-try(model.matrix(eval(frml), df), silent = T)
  if(inherits(x, "try-error")){
    if(frml[[3]] == "0"){
      # In this case, n x 0 matrix is needed
      x <- matrix(0, nrow = nrow(df), ncol = 0)
    }
  }
  if(nrow(x) < nrow(df)){ stop("Not allowed to have NAs for predictors") }

  # icenReg does not use intercepts
  if('(Intercept)' %in% colnames(x)){	
    ind = which(colnames(x) == '(Intercept)')
    x <- x[,-ind, drop = F]
  }
  ans$x <- x
  ans$xNames = colnames(x)
  
  # Getting y
  base_y = model.response(mod_frame)
  yMat = makeIntervals(base_y, mod_frame)
  ans$y = yMat
  return(ans)
}
