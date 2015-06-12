###			SEMIPARAMETRIC UTILITIES


findMaximalIntersections <- function(lower, upper){
	allVals <- sort(unique(c(lower,upper)) )
	isLeft <- allVals %in% lower
	isRight <- allVals %in% upper
	miList <- .Call("findMI", allVals, isLeft, isRight, lower, upper)
	names(miList) <- c('l_inds', 'r_inds', 'mi_l', 'mi_r')	
	return(miList)
}


fit_ICPH <- function(obsMat, covars, callText = 'ic_ph', weights){
	mi_info <- findMaximalIntersections(obsMat[,1], obsMat[,2])
	k = length(mi_info[['mi_l']])
	pmass <- rep(1/k, k)
	covars <- as.matrix(covars)
	if(callText == 'ic_ph')	fitType = as.integer(1)
	else if(callText == 'ic_po') fitType = as.integer(2)
	else stop('callText not recognized in fit_ICPH')
	myFit <- .Call('ic_sp_ch', mi_info$l_inds, mi_info$r_inds, covars, fitType, as.numeric(weights)) 
	names(myFit) <- c('p_hat', 'coefficients', 'final_llk', 'iterations', 'score')
	myFit[['T_bull_Intervals']] <- rbind(mi_info[['mi_l']], mi_info[['mi_r']])
	myFit$p_hat <- myFit$p_hat / sum(myFit$p_hat) 
	#removes occasional numerical error. Error on the order of 10^(-15), but causes problem when calculating last 
	#entry for estimates survival curve
	return(myFit)
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

getBS_coef <- function(sampDataEnv, callText = 'ic_ph'){
	xMat <- cbind(sampDataEnv$x,1)
	invertResult <- try(diag(solve(t(xMat) %*% xMat )), silent = TRUE)
	if(is(invertResult, 'try-error')) {return( rep(NA, ncol(xMat) -1) ) }
	output <- fit_ICPH(sampDataEnv$y, sampDataEnv$x, callText, sampDataEnv$w)$coefficients
	return(output)
}


expandX <- function(formula, data, fit){
	tt <- terms(fit)
	Terms <- delete.response(tt)
	 m <- model.frame(Terms, data, na.action = na.pass, xlev = fit$xlevels)
	 x <- model.matrix(Terms, m)
	 
	 return(as.matrix(x[,-1]) )
}



###		PARAMETRIC FIT UTILITIES

fit_par <- function(y_mat, x_mat, parFam = 'gamma', link = 'po', leftCen = 0, rightCen = Inf, uncenTol = 10^-6, regnames, weights){
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
	
	if(is.matrix(x_mat))			x_mat_rearranged <- rbind(x_mat[isUncen,], x_mat[isGCen,], x_mat[isLeftCen,], x_mat[isRightCen,])
	else if(length(x_mat) != 0)		x_mat_rearranged <- matrix(c(x_mat[isUncen], x_mat[isGCen], x_mat[isLeftCen], x_mat[isRightCen]), ncol = 1)
	else							x_mat_rearranged <- matrix(ncol = 0, nrow = nrow(x_mat))
	storage.mode(x_mat_rearranged) <- 'double'
	x_mat_rearranged <- as.matrix(x_mat_rearranged)	
	k_reg = ncol(x_mat_rearranged)
		
	#regnames = colnames(x_mat_rearranged)
	if(parFam == 'gamma') {parInd = as.integer(1); k_base = 2; bnames = c('log_shape', 'log_scale')}
	else if(parFam == 'weibull') {parInd = as.integer(2); k_base = 2; bnames = c('log_shape', 'log_scale')}
	else if(parFam == 'lnorm') {parInd = as.integer(3); k_base = 2; bnames = c('mu', 'log_s')}
	else if(parFam == 'exponential') {parInd = as.integer(4); k_base = 1; bnames = 'log_scale'}
	else if(parFam == 'loglogistic') {parInd = as.integer(5); k_base = 2; bnames = c('log_alpha', 'log_beta')}
	else stop('parametric family not supported')

	hessnames = c(bnames, regnames)
	
	if(link == 'po') linkInd = as.integer(1)
	else if (link == 'ph') linkInd = as.integer(2)
	else stop('link function not supported')
	
	hessian <- matrix(numeric(), nrow = (k_reg + k_base), ncol = (k_reg + k_base))
	
	fit <- .Call('ic_par', s_t, d_t, x_mat_rearranged,
				uncenInd_mat, gicInd_mat, leftCenInd, rightCenInd,
				parInd, linkInd, hessian, as.numeric(w_reordered) )
								
	names(fit) <- c('reg_pars', 'baseline', 'final_llk', 'iterations', 'hessian', 'score')
	names(fit$reg_pars) <- regnames
	names(fit$baseline) <- bnames
	colnames(fit$hessian) <- hessnames
	rownames(fit$hessian) <- hessnames
	fit$var <- -solve(fit$hessian)
	fit$coefficients <- c(fit$baseline, fit$reg_pars)
	return(fit)
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

get_etas <- function(fit, newdata = NULL){
	if(is.null(newdata)){ans <- 1; names(ans) <- 'baseline'; return(ans)}
	grpNames <- rownames(newdata)
	reducFormula <- fit$formula
	reducFormula[[2]] <- NULL
	new_x <- expandX(reducFormula, newdata, fit)
	log_etas <- as.numeric( new_x %*% fit$reg_pars)	
	etas <- exp(log_etas)
	names(etas) <- grpNames
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
get_link_fun <- function(fit){
	if(fit$model == 'po') 	return(po_link)
	if(fit$model == 'ph')	return(ph_link)
	stop('model type not recognized. Should be "ph" or "po"')
}


findUpperBound <- function(x, s_fun, link_fun, fit, eta){
	val <- 1
	fval <- 1 - link_fun(s_fun(val, fit$baseline), eta)
	tries = 0
	while(tries < 100 & fval < x){
		tries = tries + 1
		val <- val * 2
		fval <- 1 - link_fun(s_fun(val, fit$baseline), eta)
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
	fit$.dataEnv$data
}

get_tbull_mid_q <- function(p, s_t, tbulls){
	x_low <- tbulls[,1]
	x_hi <- tbulls[,2]
	x_range <- range(tbulls)
	s_t <- s_t
	k <- length(tbulls)
	
	ans <- numeric()
	
	for(i in 1:length(p)){
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
												
				if(ans[i] < 0 | ans[i] > 1) browser()								
														
			}
		}
	}
	return(ans)
}
