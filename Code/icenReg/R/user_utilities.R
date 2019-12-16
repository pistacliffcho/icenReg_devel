
vcov.icenReg_fit <- function(object,...) object$var
names.icenReg_fit <- function(x) ls(x)

#' Get Estimated Survival Curves from Semi-parametric Model for Interval Censored Data
#' 
#' @param fit model fit with \code{ic_sp} 
#' @param newdata data.frame containing covariates for which the survival curve will be fit to. 
#' Rownames from \code{newdata} will be used to name survival curve. 
#' If left blank, baseline covariates will be used
#' 
#' @description
#' Extracts the estimated survival curve(s) from an ic_sp or ic_np model for interval censored data. 
#' 
#' @details
#' Output will be a list with two elements: the first item will be \code{$Tbull_ints}, 
#' which is the Turnbull intervals. 
#' This is a k x 2 matrix, with the first column being the beginning of the 
#' Turnbull interval and the second being the end. 
#' This is necessary due to the \emph{representational non-uniqueness};
#' any survival curve that lies between the survival curves created from the
#' upper and lower limits of the Turnbull intervals will have equal likelihood. 
#' See example for proper display of this. The second item is \code{$S_curves}, 
#' or the estimated survival probability at each Turnbull interval for individuals 
#' with the covariates provided in \code{newdata}. Note that multiple rows may 
#' be provided to newdata, which will result in multiple S_curves. 
#' @author Clifford Anderson-Bergman
#' @export
getSCurves <- function(fit, newdata = NULL){
	if(inherits(fit, 'impute_par_icph'))	stop('getSCurves currently not supported for imputation model')
	if(inherits(fit, 'ic_par'))				stop('getSCurves does not support ic_par objects. Use getFitEsts() instead. See ?getFitEsts')
	etas <- get_etas(fit, newdata)
	grpNames <- names(etas)
	transFxn <- get_link_fun(fit)
	if(fit$par == 'semi-parametric' | fit$par == 'non-parametric'){
		x_l <- fit$T_bull_Intervals[1,]
		x_u <- fit$T_bull_Intervals[2,]
		x_l <- c(x_l[1], x_l)
		x_u <- c(x_l[1], x_u)
		Tbull_intervals <- cbind(x_l,  x_u)
		colnames(Tbull_intervals) <- c('lower', 'upper')
		s <- 1 - c(0, cumsum(fit$p_hat))
		ans <- list(Tbull_ints = Tbull_intervals, "S_curves" = list())
		
		for(i in 1:length(etas)){
			eta <- etas[i]
			ans[["S_curves"]][[grpNames[i] ]] <- transFxn(s, eta)
		}
		class(ans) <- 'sp_curves'
		return(ans)
	}
	else{
	  	stop('getSCurves only for semi-parametric model. Try getFitEsts')
	}
}




summary.icenReg_fit <- function(object,...)
	new('icenRegSummary', object)
summary.ic_npList <- function(object, ...)
  object

	
#' Simulates interval censored data from regression model with a Weibull baseline
#' 
#' @param n Number of samples simulated
#' @param b1 Value of first regression coefficient
#' @param b2 Value of second regression coefficient
#' @param model Type of regression model. Options are 'po' (prop. odds) and 'ph' (Cox PH)
#' @param shape shape parameter of baseline distribution
#' @param scale scale parameter of baseline distribution
#' @param inspections number of inspections times of censoring process
#' @param inspectLength max length of inspection interval
#' @param rndDigits number of digits to which the inspection time is rounded to, 
#' creating a discrete inspection time. If \code{rndDigits = NULL}, the inspection time is not rounded, 
#' resulting in a continuous inspection time
#' @param prob_cen probability event being censored. If event is uncensored, l == u
#'
#' @description
#' Simulates interval censored data from a regression model with a weibull baseline distribution. Used for demonstration
#' @details 
#' Exact event times are simulated according to regression model: covariate \code{x1} 
#' is distributed \code{rnorm(n)} and covariate \code{x2} is distributed
#' \code{1 - 2 * rbinom(n, 1, 0.5)}. Event times are then censored with a 
#' case II interval censoring mechanism with \code{inspections} different inspection times. 
#' Time between inspections is distributed as \code{runif(min = 0, max = inspectLength)}. 
#' Note that the user should be careful in simulation studies not to simulate data 
#' where nearly all the data is right censored (or more over, all the data with x2 = 1 or -1) 
#' or this can result in degenerate solutions!
#' 
#' @examples 
#' set.seed(1)
#' sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3, model = 'ph', 
#'                       shape = 2, scale = 2, inspections = 6, 
#'                       inspectLength = 1)
#' #simulates data from a cox-ph with beta weibull distribution.
#'
#' diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2, 
#'            data = sim_data, model = 'po')
#' diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2,
#'            data = sim_data, model = 'ph')
#'
#' #'ph' fit looks better than 'po'; the difference between the transformed survival
#' #function looks more constant
#' @author Clifford Anderson-Bergman
#' @export
simIC_weib <- function(n = 100, b1 = 0.5, b2 = -0.5, model = 'ph', 
					   shape = 2, scale = 2, 
					   inspections = 2, inspectLength = 2.5,
					   rndDigits = NULL, prob_cen = 1){
	rawQ <- runif(n)
    x1 <- runif(n, -1, 1)
    x2 <- 1 - 2 * rbinom(n, 1, 0.5)
    nu <- exp(x1 * b1 + x2 * b2)
    
    if(model == 'ph')		adjFun <- function(x, nu) {1 - x^(1/nu)}
	else if(model == 'po') 	adjFun <- function(x, nu) {1 - x*(1/nu) / (x * 1/nu - x + 1)}
    adjQ <- adjFun(rawQ, nu)
    trueTimes <- qweibull(adjQ, shape = shape, scale = scale)
    
    obsTimes <- runif(n = n, max = inspectLength)
    if(!is.null(rndDigits))
    	obsTimes <- round(obsTimes, rndDigits)
    
    l <- rep(0, n)
    u <- rep(0, n)
    
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    
    if(inspections > 1){
    	for(i in 2:inspections){
		    oldObsTimes <- obsTimes
    		obsTimes <- oldObsTimes + runif(n, max = inspectLength)
		    if(!is.null(rndDigits))
    			obsTimes <- round(obsTimes, rndDigits)
    		caught <- trueTimes >= oldObsTimes  & trueTimes < obsTimes
    		needsCatch <- trueTimes > obsTimes
    		u[caught] <- obsTimes[caught]
    		l[needsCatch] <- obsTimes[needsCatch]
    	}
    }
    else{
    	needsCatch <- !caught	
    }
    u[needsCatch] <- Inf
    
    if(sum(l > u) > 0)	stop('warning: l > u! Bug in code')
    
    isCensored <- rbinom(n = n, size = 1, prob = prob_cen) == 1

    l[!isCensored] <- trueTimes[!isCensored]
    u[!isCensored] <- trueTimes[!isCensored]
    
    if(sum(l == Inf) > 0){
      allTimes <- c(l,u)
      allFiniteTimes <- allTimes[allTimes < Inf]
      maxFiniteTime <- max(allFiniteTimes)
      l[l == Inf] <- maxFiniteTime
    }
    return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}


simICPO_beta <- function(n = 100, b1 = 1, b2 = -1, inspections = 1, shape1 = 2, shape2 = 2, rndDigits = NULL){
	rawQ <- runif(n)
    x1 <- rnorm(n)
    x2 <- rbinom(n, 1, 0.5) - 0.5
    nu <- exp(x1 * b1 + x2 * b2)
    adjQ <- 1 - rawQ*(1/nu) / (rawQ * 1/nu - rawQ + 1)
    trueTimes <- qbeta(adjQ, shape1 = shape1, shape2 = shape2)
    
    inspectionError = 1 / (inspections + 1)
    obsTimes <- 1 / (inspections + 1) + runif(n, min = -inspectionError, max = inspectionError)
    if(!is.null(rndDigits))
    	obsTimes <- round(obsTimes, rndDigits)
    
    l <- rep(0, n)
    u <- rep(0, n)
    
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    
    if(inspections > 1){
    	for(i in 2:inspections){
		    oldObsTimes <- obsTimes
    		obsTimes <- i / (inspections+1) + runif(n, min = -inspectionError, max = inspectionError)
		    if(!is.null(rndDigits))
    			obsTimes <- round(obsTimes, rndDigits)
    		caught <- trueTimes >= oldObsTimes  & trueTimes < obsTimes
    		needsCatch <- trueTimes > obsTimes
    		u[caught] <- obsTimes[caught]
    		l[needsCatch] <- obsTimes[needsCatch]
    	}
    }
    else{
    	needsCatch <- !caught	
    }
    u[needsCatch] <- 1
    
    if(sum(l > u) > 0)	stop('warning: l > u! Bug in code')
    
    return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}




#' Evaluate covariate effect for regression model
#' 
#' @param object      Either a formula or a model fit with \code{ic_sp} or \code{ic_par}
#' @param varName     Covariate to split data on. If left blank, will split on each covariate
#' @param data        Data. Unnecessary if \code{object} is a fit
#' @param model       Type of model. Choices are \code{'ph'} or \code{'po'} 
#' @param weights     Case weights
#' @param yType       Type of plot created. See details
#' @param factorSplit Should covariate be split as a factor (i.e. by levels)
#' @param numericCuts If \code{fractorSplit == FALSE}, cut points of covariate to stratify data on
#' @param col         Colors of each subgroup plot. If left blank, will auto pick colors
#' @param xlab        Label of x axis
#' @param ylab        Label of y axis
#' @param main        title of plot
#' @param lgdLocation Where legend should be placed. See details
#' @description Creates plots to diagnosis fit of covariate effect in a regression model. 
#' For a given variable, stratifies the data across different levels of the variable and adjusts 
#' for all the other covariates included in \code{fit} and then plots a given function to help 
#' diagnosis where covariate effect follows model assumption 
#' (i.e. either proportional hazards or proportional odds). See \code{details} for descriptions of the plots. 
#'  
#' If \code{varName} is not provided, will attempt to figure out how to divide up each covariate 
#' and plot all of them, although this may fail. 
#'@details 
#' For the Cox-PH and proportional odds models, there exists a transformation of survival curves 
#' such that the difference should be constant for subjects with different covariates. 
#' In the case of the Cox-PH, this is the log(-log(S(t|X))) transformation, for the proporitonal odds, 
#' this is the log(S(t|X) / (1 - S(t|X))) transformation. 
#'
#' The function diag_covar allows the user to easily use these transformations to diagnosis 
#' whether such a model is appropriate. In particular, it takes a single covariate and 
#' stratifies the data on different levels of that covariate. 
#' Then, it fits the semi-parametric regression model 
#' (adjusting for all other covariates in the data set) on each of these 
#' stratum and extracts the baseline survival function. If the stratified covariate does 
#' follow the regression assumption, the difference between these transformed baseline 
#' survival functions should be approximately constant. 
#'
#' To help diagnosis, the default function plotted is the transformed survival functions, 
#' with the overall means subtracted off. If the assumption holds true, then the mean 
#' removed curves should be approximately parallel lines (with stochastic noise). 
#' Other choices of \code{yType}, the function to plot, are \code{"transform"}, 
#' which is the transformed functions without the means subtracted and \code{"survival"}, 
#' which is the baseline survival distribution is plotted for each strata. 
#'
#' Currently does not support stratifying covariates that are involved in an interaction term. 
#'
#' For variables that are factors, it will create a strata for each level of the covariate, up to 20 levels. 
#' If \code{factorSplit == FALSE}, will divide up numeric covariates according to the cuts provided to numericCuts. 
#'
#' \code{lgdLocation} is an argument placed to \code{legend} dictating where the legend will be placed. 
#' If \code{lgdLocation = NULL}, will use standard placement given \code{yType}. See \code{?legend} for more details. 
#' @author Clifford Anderson-Bergman
#' @export
diag_covar <- function(object, varName, 
           data, model, weights = NULL,
           yType = 'meanRemovedTransform', 
           factorSplit = TRUE, 
           numericCuts, col, 
           xlab, ylab, main, 
           lgdLocation = NULL){
	if(!yType %in% c('survival', 'transform', 'meanRemovedTransform')) stop("yType not recognized. Options = 'survival', 'transform' or 'meanRemovedTransform'")
  if(missing(data)){
    if(!is(object, 'icenReg_fit')) stop("either object must be icenReg_fit, or formula with data supplied")
    data <- object$getRawData()
  }
  max_n_use <- nrow(data) #for backward compability. No longer need max_n_use
  subDataInfo <- subSampleData(data, max_n_use, weights)
	data <- subDataInfo$data
	weights <- subDataInfo$w
	if(is.null(weights)) weights <- rep(1, nrow(data))
	
	fullFormula <- getFormula(object)
	if(missing(model)) model <- object$model
	if(is.null(model)) stop('either object must be a fit, or model must be provided')
	if(missing(data))	data <- getData(object)
	if(missing(varName)){
		allVars <- getVarNames_fromFormula(fullFormula)
		nV <- length(allVars)
		k <- ceiling(sqrt(nV))
		if(k > 1) par(mfrow = c( ceiling(nV/k), k) )
		for(vn in allVars){
			useFactor <- length( unique((data[[vn]])) ) < 5
			diag_covar(object, vn, factorSplit = useFactor, model = model,
			           data = data, yType = yType,
			           weights = weights, lgdLocation = lgdLocation,
			           col = col)
			}
		return(invisible(NULL))
	}

	if(model == 'aft'){
	  stop('diag_covar not supported for aft model. This is because calculating
	       the non-parametric aft model is quite difficult')
	}
	
	if(model == 'ph')				s_trans <- function(x){ isOk <- x > 0 & x < 1
															ans <- numeric()
															ans[isOk] <- log(-log(x[isOk]) )
															ans
															}
	
	else if(model == 'po')		 	s_trans <- function(x){ isOk <- x > 0 & x < 1
															ans <- numeric()
															ans[isOk] <- log(x[isOk]/(1-x[isOk]))
															return(ans)
														   }
															

	allData <- data
	vals <- allData[[varName]]
	if(is.null(vals))	stop(paste('Cannot find variable', varName, 'in original dataset'))
	orgCall <- fullFormula
	newcall <- removeVarFromCall(orgCall, varName)
	if(identical(orgCall, newcall))	stop('varName not found in original call')
	
	spltInfo <- NULL
	if(factorSplit){
		factVals <- factor(vals)
		if(length(levels(factVals)) > 20) stop('Attempting to split on factor with over 20 levels. Try using numeric version of covariate and use numericCuts instead')
		spltInfo <- makeFactorSplitInfo(vals, levels(factVals))
	}
	else{
		if(missing(numericCuts))	numericCuts <- median(data[[varName]])
		spltInfo <- makeNumericSplitInfo(vals, numericCuts)
	}
	
	spltFits <- splitAndFit(newcall = newcall, data = allData, varName = varName, 
							splitInfo = spltInfo, fitFunction = ic_sp, model = model, weights = weights)
		
	allX <- numeric()
	allY <- numeric()
	fitNames <- ls(spltFits)
	for(nm in fitNames){
		allX <- c(allX, as.numeric(spltFits[[nm]]$T_bull_Intervals) )
	}
	
	xlim <- range(allX, na.rm = TRUE, finite = TRUE)
	ylim <- sort( s_trans(c(0.025, 0.975)) )
	
	if(missing(col))
		col <- 1 + 1:length(fitNames)
	if(missing(xlab))	xlab = 't'
	if(missing(main)) 	main = varName
	if(missing(ylab)){
		if(model == 'ph'){
			ylab = 'log(-log(S(t)))'
			lgdLoc <- 'bottomright'
		}
		else if(model == 'po'){
			ylab = 'log(S(t)/(1 - S(t)))'
			lgdLoc <- 'bottomleft'	
		}
		else stop('model not recognized!')
	}
	t_vals <- xlim[1] + 1:999/1000 * (xlim[2] - xlim[1])
	estList <- list()
	meanVals <- 0
	if(yType == 'survival'){
		ylab = 'S(t)'
		lgdLoc <- 'bottomleft'
		s_trans <- function(x) x
		ylim = c(0,1)
	}
	if(yType == 'meanRemovedTransform'){
		ylab = paste('Mean Removed', ylab)
	}
		
	for(i in seq_along(fitNames)){
		if(yType == 'transform' | yType == 'survival'){
			nm <- fitNames[i]
			thisCurve <- getSCurves(spltFits[[nm]])
			theseTbls <- thisCurve$Tbull_ints
			thisS <- thisCurve$S_curves$baseline
			estList[[nm]] <- list(x = theseTbls, y = s_trans(thisS) )
		}
		else if(yType == 'meanRemovedTransform'){	
			nm <- fitNames[i]
			estList[[nm]] <- s_trans(1 - getFitEsts(spltFits[[nm]], q = t_vals) )
			meanVals <- estList[[nm]] + meanVals
			}
	}

	if(yType == 'meanRemovedTransform'){
		meanVals <- meanVals/length(estList)
		ylim = c(Inf, -Inf)
		for(i in seq_along(estList)){
			theseLims <- range(estList[[i]] - meanVals, finite = TRUE, na.rm = TRUE)
			ylim[1] <- min(theseLims[1], ylim[1])
			ylim[2] <- max(theseLims[2], ylim[2])
		}
		
		yrange <- ylim[2] - ylim[1]
		ylim[2] = ylim[2] + yrange * 0.2
		ylim[1] = ylim[1] - yrange * 0.2
	}
		
	plot(NA, xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim)

	
	if(yType == 'meanRemovedTransform'){
		for(i in seq_along(estList)){
			lines(t_vals, estList[[i]] - meanVals, col = col[i])
		}
	}	
	else if(yType == 'transform' | yType == 'survival'){
		for(i in seq_along(estList)){
			xs <- estList[[i]]$x
			y  <- estList[[i]]$y
			lines(xs[,1], y, col = col[i], type = 's')
			lines(xs[,2], y, col = col[i], type = 's')
		}
	}
	if(!is.null(lgdLocation))	lgdLoc <- lgdLocation
	legend(lgdLoc, legend = fitNames, col = col, lwd = 1)
}




#' Get Survival Curve Estimates from icenReg Model
#' 
#' @param fit      model fit with \code{ic_par} or \code{ic_sp}
#' @param newdata  \code{data.frame} containing covariates
#' @param p        Percentiles
#' @param q        Quantiles
#' @description 
#' Gets probability or quantile estimates from a \code{ic_sp}, \code{ic_par} or \code{ic_bayes} object. 
#' Provided estimates conditional on regression parameters found in \code{newdata}. 
#' @details
#' For the \code{ic_sp} and \code{ic_par}, the MLE estimate is returned. For \code{ic_bayes}, 
#' the MAP estimate is returned. To compute the posterior means, use \code{sampleSurv}.
#' 
#' If \code{newdata} is left blank, baseline estimates will be returned (i.e. all covariates = 0). 
#' If \code{p} is provided, will return the estimated F^{-1}(p | x). If \code{q} is provided, 
#' will return the estimated F(q | x). If neither \code{p} nor \code{q} are provided, 
#' the estimated conditional median is returned.
#'  
#' In the case of \code{ic_sp}, the MLE of the baseline survival is not necessarily unique, 
#' as probability mass is assigned to disjoint Turnbull intervals, but the likelihood function is 
#' indifferent to how probability mass is assigned within these intervals. In order to have a well 
#' defined estimate returned, we assume probability is assigned uniformly in these intervals. 
#' In otherwords, we return *a* maximum likelihood estimate, but don't attempt to characterize *all* maximum 
#' likelihood estimates with this function. If that is desired, all the information needed can be 
#' extracted with \code{getSCurves}.
#' @examples 
#' simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
#' inspections = 6, inspectLength = 1)
#' fit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2,
#'              data = simdata)
#' new_data <- data.frame(x1 = c(1,2), x2 = c(-1,1))
#' rownames(new_data) <- c('grp1', 'grp2')
#' 
#' estQ <- getFitEsts(fit, new_data, p = c(.25, .75))
#' 
#' estP <- getFitEsts(fit, q = 400)
#' @author Clifford Anderson-Bergman
#' @export
getFitEsts <- function(fit, newdata = NULL, p, q){
  if(is.null(newdata)){}
  else{
    if(!identical(newdata, "midValues")){ 
      if(!is(newdata, "data.frame")){ 
        stop("newdata should be a data.frame")
      }
    }
  }
  etas <- get_etas(fit, newdata)
  
  
  if(missing(p))	p <- NULL
  if(missing(q))  q <- NULL
  if(!is.null(q)) {xs <- q; type = 'q'}
  else{ 
    type = 'p'
    if(is.null(p)) xs <- 0.5
    else		   xs <- p
  }
  
  if(length(etas) == 1){etas <- rep(etas, length(xs))}
  if(length(xs) == 1){xs <- rep(xs, length(etas))}
  if(length(etas) != length(xs) ) stop('length of p or q must match nrow(newdata) OR be 1')

  regMod <- fit$model
  
  if(inherits(fit, 'sp_fit'))	{
    scurves <- getSCurves(fit, newdata = NULL)
    baselineInfo <- list(tb_ints = scurves$Tbull_ints, s = scurves$S_curves$baseline)
    baseMod = 'sp'
  }
  if(inherits(fit, 'par_fit') | inherits(fit, 'bayes_fit')){	
    baseMod <- fit$par
    baselineInfo <- fit$baseline
  }
  
  if(fit$model == 'po' | fit$model == 'ph' | fit$model == 'none'){
    surv_etas <- etas
    scale_etas <- rep(1, length(etas)) 
  }
  
  else if(fit$model == 'aft'){
    scale_etas <- etas
    surv_etas <- rep(1, length(etas)) 
  }
  else stop('model not recognized in getFitEsts')
  
  if(type == 'q'){
    ans <- getSurvProbs(xs /scale_etas, surv_etas, 
                        baselineInfo = baselineInfo, 
                        regMod = regMod, baseMod = baseMod)
#    ans <- ans * scale_etas
    return(ans)
  }
  else if(type == 'p'){
#    xs <- xs / scale_etas
    ans <- getSurvTimes(xs, surv_etas, 
                        baselineInfo = baselineInfo, 
                        regMod = regMod, baseMod = baseMod)
    ans <- ans * scale_etas
  return(ans)
  }
}


#' Compare parametric baseline distributions with semi-parametric baseline
#' 
#' @param object       Either a formula or a model fit with \code{ic_sp} or \code{ic_par}
#' @param data         Data. Unnecessary if \code{object} is a fit
#' @param model        Type of model. Choices are \code{'ph'} or \code{'po'}
#' @param dists        Parametric baseline fits	
#' @param cols         Colors of baseline distributions
#' @param weights      Case weights
#' @param lgdLocation  Where legend will be placed. See \code{?legend} for more details
#' @param useMidCovars Should the distribution plotted be for covariates = mean values instead of 0
#'
#' @description 
#' Creates plots to diagnosis fit of different choices of parametric baseline model. 
#' Plots the semi paramtric model against different choices of parametric models. 
#' 
#' @details
#' If \code{useMidCovars = T}, then the survival curves plotted are for fits with the mean covariate value, 
#' rather than 0. This is because often the baseline distribution (i.e. with all covariates = 0) will be 
#' far away from the majority of the data.
#' 
#' @examples data(IR_diabetes)
#' fit <- ic_par(cbind(left, right) ~ gender, 
#'              data = IR_diabetes)
#'
#' diag_baseline(fit, lgdLocation = "topright", 
#'              dist = c("exponential", "weibull", "loglogistic"))
#'
#' @author Clifford Anderson-Bergman
#' @export
diag_baseline <- function(object, data, model = 'ph', weights = NULL,
						  dists = c('exponential', 'weibull', 'gamma', 'lnorm', 'loglogistic', 'generalgamma'),
						  cols = NULL, lgdLocation = 'bottomleft',
						  useMidCovars = T){
  if(model == 'aft'){
    stop('diag_baseline not supported for aft model. This is because calculating
	       the non-parametric aft model is quite difficult')
  }  
	newdata = NULL
	if(useMidCovars) newdata <- 'midValues'
	formula <- getFormula(object)
	if(missing(data))	data <- getData(object)
	max_n_use = nrow(data)	#no longer necessary, for backward compatability				  	
	
	subDataInfo <- subSampleData(data, max_n_use, weights)
	sp_data <- subDataInfo$data
	weights <- subDataInfo$w

	sp_fit <- ic_sp(formula, data = sp_data, bs_samples = 0, model = model)
	plot(sp_fit, newdata)
	xrange <- range(getSCurves(sp_fit)$Tbull_ints, finite = TRUE)
	grid <- xrange[1] + 0:100/100 *(xrange[2] - xrange[1])
	if(is.null(cols)) cols <- 1 + 1:length(dists)
	for(i in seq_along(dists)){
		this_dist <- dists[i]
		par_fit <- ic_par(formula, data = data, model = model, dist = this_dist)
		y <- getFitEsts(par_fit, newdata = newdata, q = grid)
		lines(grid, 1 - y, col = cols[i])
	}
	legend(lgdLocation, legend = c('Semi-parametric', dists), col = c('black', cols), lwd = 1)
}

#' Predictions from icenReg Regression Model
#' 
#' @param object   Model fit with \code{ic_par} or \code{ic_sp}
#' @param type     type of prediction. Options include \code{"lp", "response"}
#' @param newdata  \code{data.frame} containing covariates
#' @param ...      other arguments (will be ignored)
#'
#' @description   
#' Gets various estimates from an \code{ic_np}, \code{ic_sp} or \code{ic_par} object.
#' @details 
#' If \code{newdata} is left blank, will provide estimates for original data set. 
#'
#' For the argument \code{type}, there are two options. \code{"lp"} provides the 
#' linear predictor for each subject (i.e. in a proportional hazards model, 
#' this is the log-hazards ratio, in proportional odds, the log proporitonal odds), 
#' \code{"response"} provides the median response value for each subject, 
#' *conditional on that subject's covariates, but ignoring their actual response interval*. 
#' Use \code{imputeCens} to impute the censored values.
#' @examples 
#' simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
#'                       inspections = 6,
#'                       inspectLength = 1)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' imputedValues <- predict(fit)
#' @author Clifford Anderson-Bergman
#' @export
predict.icenReg_fit <- function(object, type = 'response',
                                newdata = NULL, ...)
      #imputeOptions = fullSample, fixedParSample, median
  {
  if(is.null(newdata)) newdata <- object$getRawData()
  if(type == 'lp')
    return( log(get_etas(object, newdata = newdata)))
  if(type == 'response')
    return(getFitEsts(fit = object, 
                      newdata = newdata, 
                      p = 0.5))
  stop('"type" not recognized: options are "lp", "response" and "impute"')
}

#' Impute Interval Censored Data from icenReg Regression Model
#' 
#' @param fit         icenReg model fit 
#' @param newdata     \code{data.frame} containing covariates and censored intervals. If blank, will use data from model
#' @param imputeType  type of imputation. See details for options
#' @param samples  Number of imputations (ignored if \code{imputeType = "median"}) 
#' 
#' @description
#' Imputes censored responses from data. 
#' @details 	
#'  If \code{newdata} is left blank, will provide estimates for original data set. 
#' 
#'  There are several options for how to impute. \code{imputeType = 'median'} 
#'  imputes the median time, conditional on the response interval, covariates and 
#'  regression parameters at the MLE. To get random imputations without accounting
#'  for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a 
#'  random sample of the response variable, conditional on the response interval, 
#'  covariates and estimated parameters at the MLE. Finally, 
#'  \code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
#'  (assuming asymptotic normality) and then takes a random sample 
#'  of the response variable, conditional on the response interval, 
#'  covariates, and the random sample of the coefficients. 
#'  
#' @examples 
#' simdata <- simIC_weib(n = 500)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' imputedValues <- imputeCens(fit)
#' @author Clifford Anderson-Bergman
#' @export
imputeCens<- function(fit, newdata = NULL, imputeType = 'fullSample', samples = 5){
  if(is.null(newdata)) newdata <- fit$getRawData()
  yMat <- expandY(fit$formula, newdata, fit)
  p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
  p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
  ans <- matrix(nrow = length(p1), ncol = samples)
  storage.mode(ans) <- 'double'
  if(imputeType == 'median'){
    isBayes <- is(fit, 'bayes_fit')
    if(isBayes){
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    p_med <- (p1 + p2)/2
    ans <- getFitEsts(fit, newdata, p = p_med)
    isLow <- ans < yMat[,1]
    ans[isLow] <- yMat[isLow,1]
    isHi <- ans > yMat[,2]
    ans[isHi] <- yMat[isHi]
    if(isBayes) setSamplablePars(fit, orgCoefs)
#    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if(imputeType == 'fixedParSample'){
    isBayes <- is(fit, 'bayes_fit')
    if(isBayes){
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    for(i in 1:samples){
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      rownames(ans) <- rownames(newdata)
      ans <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    if(isBayes) setSamplablePars(fit, orgCoefs)
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if(imputeType == 'fullSample'){
    isSP <- is(fit, 'sp_fit')
    isBayes <- is(fit, 'bayes_fit')
    for(i in 1:samples){
      orgCoefs <- getSamplablePars(fit)
      if(isBayes){
        sampledCoefs <- sampBayesPar(fit)
      }
      else if(!isSP){
        coefVar <- getSamplableVar(fit)
        sampledCoefs <- sampPars(orgCoefs, coefVar)
      }
      else{
        sampledCoefs <- getBSParSample(fit)
      }
      setSamplablePars(fit, sampledCoefs)
      p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
      p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      fastMatrixInsert(theseImputes, ans, colNum = i)
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  stop('imputeType type not recognized.')
}

#' Draw samples from an icenReg model
#' 
#' @param fit         icenReg model fit 
#' @param newdata     \code{data.frame} containing covariates. If blank, will use data from model
#' @param sampleType  type of samples See details for options
#' @param samples     Number of samples 
#' 
#' @description
#' Samples response values from an icenReg fit conditional on covariates, but not 
#' censoring intervals. To draw response samples conditional on covariates and 
#' restrained to intervals, see \code{imputeCens}.
#'  
#'  
#' @details 	
#'  Returns a matrix of samples. Each row of the matrix corresponds with a subject with the 
#'  covariates of the corresponding row of \code{newdata}. For each column of the matrix, 
#'  the same sampled parameters are used to sample response variables. 
#'  
#'  If \code{newdata} is left blank, will provide estimates for original data set. 
#' 
#'  There are several options for how to sample. To get random samples without accounting
#'  for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a 
#'  random sample of the response variable, conditional on the response interval, 
#'  covariates and estimated parameters at the MLE. Alternatively, 
#'  \code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
#'  (assuming asymptotic normality for the ic_par) and then takes a random sample 
#'  of the response variable, conditional on the response interval, 
#'  covariates, and the random sample of the coefficients. 
#'  
#' @examples 
#' simdata <- simIC_weib(n = 500)
#'
#' fit <- ic_par(cbind(l, u) ~ x1 + x2,
#'               data = simdata)
#'
#' newdata = data.frame(x1 = c(0, 1), x2 = c(1,1))
#'
#' sampleResponses <- ic_sample(fit, newdata = newdata, samples = 100)
#' @author Clifford Anderson-Bergman
#' @export
ic_sample <- function(fit, newdata = NULL, sampleType = 'fullSample',
                      samples = 5){
  if(is.null(newdata)) newdata <- fit$getRawData()
  yMat <- cbind(rep(-Inf, nrow(newdata)), rep(Inf, nrow(newdata)))
#  p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
#  p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
  nRow = nrow(newdata)
  p1 = rep(0, nrow(newdata))
  p2 = rep(1, nrow(newdata))
  ans <- matrix(nrow = length(p1), ncol = samples)
  storage.mode(ans) <- 'double'
  isSP <- is(fit, 'sp_fit')
  isBayes <- is(fit, 'bayes_fit')
  if(sampleType == 'fixedParSample'){
    if(isBayes){
      orgCoefs <- getSamplablePars(fit)
      map_ests <- c(fit$MAP_baseline, fit$MAP_reg_pars)
      setSamplablePars(fit, map_ests)
    }
    for(i in 1:samples){
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      ans[,i] <- fastMatrixInsert(theseImputes, ans, colNum = i)
    }
    if(isBayes) setSamplablePars(fit, orgCoefs)
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  if(sampleType == 'fullSample'){
    for(i in 1:samples){
      orgCoefs <- getSamplablePars(fit)
      if(isBayes){ sampledCoefs <- sampBayesPar(fit) }
      else if(!isSP){
        coefVar <- getSamplableVar(fit)
        sampledCoefs <- sampPars(orgCoefs, coefVar)
      }
      else{ sampledCoefs <- getBSParSample(fit) }
      setSamplablePars(fit, sampledCoefs)
#      p1 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,1]) ) 
#      p2 <- getFitEsts(fit, newdata, q = as.numeric(yMat[,2]) ) 
      p_samp <- runif(length(p1), p1, p2)
      theseImputes <- getFitEsts(fit, newdata, p = p_samp)
      isLow <- theseImputes < yMat[,1]
      theseImputes[isLow] <- yMat[isLow,1]
      isHi <- theseImputes > yMat[,2]
      theseImputes[isHi] <- yMat[isHi,2]
      ans[,i] <- theseImputes
#      fastMatrixInsert(theseImputes, ans, colNum = i)
      setSamplablePars(fit, orgCoefs)
    }
    rownames(ans) <- rownames(newdata)
    return(ans)
  }
  stop('sampleType type not recognized.')
}


sampleSurv_slow <- function(fit, newdata, p = NULL, q = NULL, samples = 100){
  if(nrow(newdata) > 1) stop('newdata must be a single row')
  # Checking what type of look up to do 
  input_type = NULL
  if(!is.null(p)){
    input_type = 'p'
    nCol = length(p)
  }
  if(!is.null(q)){
    input_type = 'q'
    nCol = length(q)
  }
  if(is.null(input_type)){ 
    input_type = 'p'
    p = c(0.1, .25, .5, .75, .9)
    nCol = length(p)
  }
  is_np <- is(fit, 'sp_fit') | is(fit, 'np_fit')
  if(is_np) stop("Can only sample survival estimates for parametric/Bayesian models")
  isBayes <- is(fit, 'bayes_fit')
  ans <- matrix(nrow = samples, ncol = nCol)
  for(i in 1:samples){
    orgCoefs <- getSamplablePars(fit)
    if(isBayes){ sampledCoefs <- sampBayesPar(fit) }
    else{
      coefVar <- getSamplableVar(fit)
      sampledCoefs <- sampPars(orgCoefs, coefVar)
    }
    setSamplablePars(fit, sampledCoefs)
    if(input_type == 'p') theseSamples <- getFitEsts(fit, newdata, p = p)
    else theseSamples <- getFitEsts(fit, newdata, q = q)
    ans[i,] <- theseSamples
    setSamplablePars(fit, orgCoefs)
  }
  if(input_type == 'p') these_colnames <- paste("p =", p)
  else these_colnames <- paste("t =", q)
  return(ans)
}

#' Samples fitted survival function
#' 
#' @param fit       Either an ic_bayes or ic_par fit
#' @param newdata   A data.frame with a single row of covariates
#' @param p         A set of survival probabilities to sample corresponding time for
#' @param q         A set of times to sample corresponding cumulative probability for
#' @param samples   Number of samples to draw
#' @details For Bayesian models, draws samples from the survival distribution with a given set of covariates.
#' Does this by first drawing a set of parameters (both regression and baseline) from \code{fit$samples} and then computing the quantiles of 
#' the distribution (if \code{p} is provided) or the CDF at \code{q}. 
#' 
#' If a \code{ic_par} model is provided, the procedure is the same, but the sampled parameters are drawn using
#' the normal approximation.
#' 
#' Not compatible with \code{ic_np} or \code{ic_sp} objects.
#' @author Clifford Anderson-Bergman
#' 
#' @examples 
#' data("IR_diabetes")
#' fit <- ic_par(cbind(left, right) ~ gender, data = IR_diabetes)
#' 
#' newdata <- data.frame(gender = "male")
#' time_samps <- sampleSurv(fit, newdata, 
#'                          p = c(0.5, .9), 
#'                          samples = 100)
#' # 100 samples of the median and 90th percentile for males                        
#' 
#' prob_samps <- sampleSurv(fit, newdata, 
#'                          q = c(10, 20),
#'                          samples = 100)
#' # 100 samples of the cumulative probability at t = 10 and 20 for males                        
#' @export
sampleSurv <- function(fit, newdata = NULL, p = NULL, q = NULL, samples = 100){
  if(!is.null(newdata)){
    if(nrow(newdata) > 1) stop('newdata must be a single row')
  }
  # Checking what type of look up to do 
  input_type = NULL
  if(!is.null(p)){
    input_type = 'p'
    input = p
    nCol = length(p)
  }
  if(!is.null(q)){
    input_type = 'q'
    input = q
    nCol = length(q)
  }
  if(is.null(input_type)){ 
    input_type = 'p'
    input = c(0.1, .25, .5, .75, .9)
    nCol = length(input)
  }
  is_np <- is(fit, 'sp_fit') | is(fit, 'np_fit')
  if(is_np) stop("Can only sample survival estimates for parametric/Bayesian models")
  isBayes <- is(fit, 'bayes_fit')
  ans <- matrix(nrow = samples, ncol = nCol)
  etas_and_base <- sample_etas_and_base(fit, samples = samples, newdata = newdata)
  etas <- etas_and_base$etas
  if(length(etas) == 1) etas <- rep(etas, samples)
  baseMat <- etas_and_base$baseMat
  if(nrow(baseMat) != samples) stop("nrow(baseMat) != samples")
  if(length(etas)  != samples) stop("length(etas) != samples")
  for(i in 1:nCol){
    this_input = rep(input[i], samples)
    if(input_type == 'p'){ 
      ans[,i] <- computeConditional_q(this_input, 
                                      etas,
                                      baseMat, 
                                      fit$model,
                                      fit$par) 
    }
    if(input_type == 'q'){ 
      ans[,i] <- computeConditional_p(this_input, 
                                      etas,
                                      baseMat, 
                                      fit$model,
                                      fit$par) 
    }
    
  }
  if(input_type == 'p') these_colnames <- paste("Q(p) =", input)
  else these_colnames <- paste("F(t) =", input)
  colnames(ans) <- these_colnames
  return(ans)
}




dGeneralGamma <- function(x, mu, s, Q){
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('dGeneralGamma', x, mu, s, Q)
  return(ans)
}

qGeneralGamma <- function(p, mu, s, Q){
  x <- p
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('qGeneralGamma', x, mu, s, Q)
  return(ans)
  
}

pGeneralGamma <- function(q, mu, s, Q){
  x <- q
  max_n <- getMaxLength(list(x, mu, s, Q) )
  x <- updateDistPars(x, max_n)
  mu <- updateDistPars(mu, max_n)
  s <- updateDistPars(s, max_n)
  Q <- updateDistPars(Q, max_n)
  
  ans <- .Call('qGeneralGamma', x, mu, s, Q)
  return(ans)
  
}


icqqplot <- function(par_fit){
 spfit <- makeQQFit(par_fit)
 baseSCurves <- getSCurves(spfit)
 baseS <- baseSCurves$S_curves$baseline
 baseUpper <- baseSCurves$Tbull_ints[,2]
 baseLower <- baseSCurves$Tbull_ints[,1]
 
 parUpperS <- 1 - getFitEsts(fit = par_fit, q = baseUpper)
 parLowerS <- 1 - getFitEsts(fit = par_fit, q = baseLower)
 
 plot(baseS, parUpperS, xlim = c(0,1), ylim = c(0,1), main = 'QQ Plot', xlab = c('Unconstrained Baseline Survival'),
      ylab = 'Parametric Survival', type = 's')
 lines(baseS, parLowerS, type = 's')
 lines(c(0,1), c(0,1), col = 'red')
}


#' Convert current status data into interval censored format
#' 
#' @param time           Time of inspection
#' @param eventOccurred  Indicator if event has occurred. 0/FALSE implies event has not occurred by inspection time. 
#' 
#' @details Converts current status data to the interval censored format for 
#' usage in icenReg.
#' 
#' @examples
#' 
#' simData <- simCS_weib()
#' # Simulate current status data
#' 
#' head(cs2ic(simData$time, simData$event))
#' # Converting data to current status format
#' 
#' fit <- ic_par(cs2ic(time, event) ~ x1 + x2, data = simData)
#' # Can be used directly in formula
#' 
#' @export
cs2ic <- function(time, 
                  eventOccurred){
  if(!all(eventOccurred == 0 | eventOccurred == 1)){stop("eventOccurred must be all 1's and 0's or TRUE's and FALSE's")} 
  if(length(time) != length(eventOccurred)) stop('length(time) != length(eventOccurred)')
  l = rep(0, length(time))
  u = rep(Inf, length(time))
  logical_event <- eventOccurred == 1
  l[!eventOccurred] <- time[!eventOccurred]
  u[eventOccurred]  <- time[eventOccurred]
  return(cbind(l, u))
}

#' Confidence/Credible intervals for survival curves
#' 
#' @param fit Fitted model from \code{ic_par} or \code{ic_bayes}
#' @param p Percentiles of distribution to sample
#' @param q Times of disitribution to sample. Only p OR q should be specified, not p AND q
#' @param newdata \code{data.frame} containing covariates for survival curves
#' @param ci_level Confidence/credible level
#' @param MC_samps Number of Monte Carlo samples taken
#' @details Creates a set of confidence intervals for the survival curves conditional
#' on the covariates provided in \code{newdata}. Several rows can be provided in newdata;
#' this will lead to several sets of confidence/credible intervals. 
#' 
#' For Bayesian models, these are draws directly from the posterior; a set of parameters drawn 
#' from those saved in \code{fit$samples} repeatedly and then for each set of parameters, 
#' the given set of quantiles is calculated. For parametric models, the procedure is virtually the 
#' same, but rather than randomly drawing rows from saved samples, random samples are drawn using
#' the asymptotic normal approximation of the estimator. 
#' 
#' This function is not compatible with \code{ic_np} or \code{ic_sp} objects, as the distribution 
#' of the baseline distribution of these estimators is still an open question. 
#' @author Clifford Anderson-Bergman
#' @examples 
#' data("IR_diabetes")
#' fit <- ic_par(cbind(left, right) ~ gender, 
#'                 data = IR_diabetes)
#' 
#' # Getting confidence intervals for survival curves
#' # for males and females
#' newdata <- data.frame(gender = c("male", "female"))
#' rownames(newdata) <- c("Males", "Females")
#' diab_cis <- survCIs(fit, newdata)
#' diab_cis
#' 
#' # Can add this to any plot
#' plot(fit, newdata = newdata, 
#'      cis = FALSE)
#' # Would have been included by default
#' lines(diab_cis, col = c("black", "red"))
#' @export
survCIs <- function(fit, newdata = NULL, 
                    p = NULL, 
                    q = NULL, 
                    ci_level = 0.95,
                    MC_samps = 4000){
  if(is.null(p) & is.null(q)){ p = c(0:19 * .05 + 0.025) }
  if(!is.null(p) & !is.null(q)){ stop('Need to provide p OR q, not both') }
  if(!(inherits(fit, 'par_fit') | inherits(fit, 'bayes_fit')))
     stop("fit must be either from ic_par() or ic_bayes()")
  ans <- surv_cis(fit = fit, newdata = newdata, 
                  p = p, q = q, ci_level = ci_level, 
                  MC_samps = MC_samps)
  return(ans)
}

