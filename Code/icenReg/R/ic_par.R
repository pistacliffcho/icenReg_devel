#' Parametric Regression  Models for Interval Censored Data
#' 
#' @param formula     Regression formula. Response must be a \code{Surv} object of type
#'  \code{'interval2'} or \code{cbind}. See details.
#' @param data        Dataset
#' @param model       What type of model to fit. Current choices are "\code{ph}" (proportional hazards), 
#' "\code{po}" (proportional odds) or "\code{aft}" (accelerated failure time)
#' @param dist        What baseline parametric distribution to use. See details for current choices
#' @param weights     vector of case weights. Not standardized; see details
#'
#' @description Fits a parametric regression model for interval censored data. 
#' Can fita proportional hazards, proportional odds or accelerated failure time model.  
#'
#' @details Currently supported distributions choices are "exponential", "weibull", "gamma", 
#' "lnorm", "loglogistic" and "generalgamma" (i.e. generalized gamma distribution). 
#'
#' Response variable should either be of the form \code{cbind(l, u)} or \code{Surv(l, u, type = 'interval2')}, 
#' where \code{l} and \code{u} are the lower and upper ends of the interval known to contain the event of interest. 
#' Uncensored data can be included by setting \code{l == u}, right censored data can be included by setting 
#' \code{u == Inf} or \code{u == NA} and left censored data can be included by setting \code{l == 0}.
#'
#' Does not allow uncensored data points at t = 0 (i.e. \code{l == u == 0}), as this will 
#' lead to a degenerate estimator for most parametric families. Unlike the current implementation 
#' of survival's \code{survreg}, does allow left side of intervals of positive length to 0 and 
#' right side to be \code{Inf}. 
#'
#' In regards to weights, they are not standardized. This means that if weight[i] = 2, 
#' this is the equivalent to having two observations with the same values as subject i. 
#' 
#' 
#' For numeric stability, if abs(right - left) < 10^-6, observation are considered 
#' uncensored rather than interval censored with an extremely small interval. 
#' @examples
#' data(miceData)
#'
#' logist_ph_fit <- ic_par(Surv(l, u, type = 'interval2') ~ grp, 
#'                        data = miceData, dist = 'loglogistic')
#' 
#' logist_po_fit <- ic_par(cbind(l, u) ~ grp, 
#'                         data = miceData, dist = 'loglogistic',
#'                        model = 'po')
#'
#' summary(logist_ph_fit)
#' summary(logist_po_fit)
#' @author Clifford Anderson-Bergman
#' @export
ic_par <- function(formula, data, model = 'ph', dist = 'weibull', weights = NULL){
  if(missing(data)) data <- environment(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  #    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts)
  if(is.matrix(x))	xNames <- colnames(x)
  else				xNames <- as.character(formula[[3]])
  if('(Intercept)' %in% colnames(x)){	
    ind = which(colnames(x) == '(Intercept)')
    x <- x[,-ind]
    xNames <- xNames[-ind]
  }
  
  yMat <- as.matrix(y)[,1:2]
  
  if(is(y, "Surv")){
    rightCens <- mf[,1][,3] == 0
    yMat[rightCens,2] <- Inf
    
    exact <- mf[,1][,3] == 1
    yMat[exact, 2] = yMat[exact, 1]
  }
  storage.mode(yMat) <- 'double'
  
  if(sum(is.na(mf)) > 0)
    stop("NA's not allowed. If this is supposed to be right censored (i.e. [4, NA] was supposed to be right censored at t = 4), replace NA with Inf")
  
  testMat <- cbind(x, 1)
  invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
  if(is(invertResult, 'try-error'))
    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
  
  callText <- paste(dist, model)
  
  if(is.null(weights))	weights = rep(1, nrow(yMat))
  if(length(weights) != nrow(yMat))	stop('weights improper length!')
  if(min(weights) < 0)				stop('negative weights not allowed!')
  if(sum(is.na(weights)) > 0)			stop('cannot have weights = NA')
  if(is.null(ncol(x))) recenterCovar = FALSE
  fitInfo <- fit_par(yMat, x, parFam = dist, link = model, 
                     leftCen = 0, rightCen = Inf, uncenTol = 10^-6, 
                     regnames = xNames, weights = weights,
                     callText = callText)
  fitInfo$call = cl
  fitInfo$formula = formula
  fitInfo$.dataEnv = new.env()
  if(!missing(data)){ fitInfo$.dataEnv$data = data }
  fitInfo$par = dist
  fitInfo$model = model
  fitInfo$terms <- mt
  fitInfo$xlevels <- .getXlevels(mt, mf)
  
  return(fitInfo)
}



fit_par <- function(y_mat, x_mat, parFam = 'gamma', link = 'po', 
                    leftCen = 0, rightCen = Inf, 
                    uncenTol = 10^-6, regnames, 
                    weights, callText){
  
  parList<- make_par_fitList(y_mat, x_mat, parFam = parFam, 
                             link = link, leftCen = leftCen, rightCen = rightCen,
                             uncenTol = uncenTol, regnames = regnames,
                             weights = weights, callText = callText)
  
  c_fit <- ic_parList(parList)
  
  
  fit <- new(callText)
  fit$reg_pars      <- c_fit$reg_pars
  fit$baseline      <- c_fit$baseline
  fit$llk           <- c_fit$llk
  fit$iterations    <- c_fit$iterations
  fit$hessian       <- c_fit$hessian
  fit$score         <- c_fit$score
  fit$par           <- parFam
  
  recenterCovar <- FALSE
#   if(recenterCovar == TRUE){
#     fit$pca_coefs <- fit$reg_pars
#     fit$pca_hessian  <- fit$hessian
#     fit$pca_info <- prcomp_xmat
#     
#     allPars <- c(fit$baseline, fit$reg_pars)
#     
#     transformedPar <- PCAFit2OrgParFit(prcomp_xmat, fit$pca_hessian, allPars, k_base)
#     fit$baseline   <- transformedPar$pars[1:k_base]
#     fit$reg_pars   <- transformedPar$pars[-1:-k_base]
#     fit$var        <- transformedPar$var	
#     fit$hessian    <- solve(fit$var)
#     fit$baseOffset <- as.numeric(fit$reg_pars %*% prcomp_xmat$center)
#   }
  
  names(fit$reg_pars)    <- parList$regnames
  names(fit$baseline)    <- parList$bnames
  colnames(fit$hessian)  <- parList$hessnames
  rownames(fit$hessian)  <- parList$hessnames
  if(recenterCovar == FALSE){
    fit$var <- -solve(fit$hessian)
    fit$baseOffset = 0
  }
  fit$coefficients <- c(fit$baseline, fit$reg_pars)
  return(fit)
}

