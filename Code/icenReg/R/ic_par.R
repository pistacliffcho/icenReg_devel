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
  checkFor_cluster(formula)
  
  # Information about orginal call to function. Useful for expanding X in predict(fit, newdata)
  call_base = match.call(expand.dots = FALSE)
  call_info = readingCall(call_base)
  
  reg_items = make_xy(formula, data)
  yMat <- reg_items$y
  x <- reg_items$x
  xNames = reg_items$xNames

  testMat <- cbind(x, 1)
  invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
  if(is(invertResult, 'try-error'))
    stop('covariate matrix is computationally singular! Make sure not to add intercept to model, also make sure every factor has observations at every level')
  
  callText <- paste(dist, model)
  
  if(is.null(weights))	weights = rep(1, nrow(yMat))
  if(length(weights) != nrow(yMat))	stop('weights improper length!')
  if(min(weights) < 0)				stop('negative weights not allowed!')
  if(sum(is.na(weights)) > 0)			stop('cannot have weights = NA')

  # Recentering covariates
  covarOffset <- icColMeans(x)
  x <- t(t(x) - covarOffset)
  
  fitInfo <- fit_par(yMat, x, parFam = dist, link = model, 
                     leftCen = 0, rightCen = Inf, uncenTol = 10^-6, 
                     regnames = xNames, weights = weights,
                     callText = callText)
  fitInfo$call = call_base
  fitInfo$formula = formula
  fitInfo$.dataEnv = new.env()
  if(!missing(data)){ fitInfo$.dataEnv$data = data }
  fitInfo$par = dist
  fitInfo$model = model
  fitInfo$terms <- call_info$mt
  fitInfo$xlevels <- .getXlevels(call_info$mt, call_info$mf)
  fitInfo$covarOffset <- matrix(covarOffset, nrow = 1)
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
  

  names(fit$reg_pars)    <- parList$regnames
  names(fit$baseline)    <- parList$bnames
  colnames(fit$hessian)  <- parList$hessnames
  rownames(fit$hessian)  <- parList$hessnames
  fit$var                <- -solve(fit$hessian)
  fit$coefficients <- c(fit$baseline, fit$reg_pars)
  return(fit)
}

icColMeans <- function(x){
  dims <- dim(x)
  if(is.null(dims)) return(mean(x))
  if(length(dims) == 2) return(colMeans(x))
  stop('data type unrecognized')
}