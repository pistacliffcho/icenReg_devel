#' Bayesian Regression  Models for Interval Censored Data
#' 
#' @param formula        Regression formula. Response must be a \code{Surv} object of type
#'  \code{'interval2'} or \code{cbind}. See details.
#' @param data           Dataset
#' @param model          What type of model to fit. Current choices are "\code{ph}" (proportional hazards), 
#' "\code{po}" (proportional odds) or "\code{aft}" (accelerated failure time)
#' @param dist           What baseline parametric distribution to use. See details for current choices
#' @param weights        vector of case weights. Not standardized; see details
#' @param logPriorFxn    An R function that computes the log prior
#' @param controls       Control parameters passed to samplers 
#' @param useMCores      Should multiple cores be used? Each core is used to run a single chain.
#' 
#' @description Fits a Bayesian regression model for interval censored data. 
#' Can fit a proportional hazards, proportional odds or accelerated failure time model.  
#'
#' @details Currently supported distributions choices are "exponential", "weibull", "gamma", 
#' "lnorm", "loglogistic" and "generalgamma" (i.e. generalized gamma distribution). 
#'
#' The \code{logPriorFxn} should take in the a vector of values corresponding to \emph{all}
#' the parameters of the model (baseline parameters first, regression parameters second) and returns the
#' log prior, calculated up to an additive constant. Default behavior is to use a flat prior. 
#' See examples for an example of using the log prior function.
#'
#' Sampling is done by a single MH block updater on all the parameters. 
#' See \code{?bayesControls} for more details. 
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
#' flat_prior_model <- ic_bayes(cbind(l, u) ~ grp, data = miceData)
#' # Default behavior is flat prior
#' 
#' priorFxn <- function(pars){
#'  ans <- 0
#'  ans <- ans + dnorm(pars[1], log = TRUE)
#'  ans <- ans + dnorm(pars[3], sd = 0.25, log = TRUE)
#' }
#' # Prior function puts N(0,1) prior on baseline shape parameter (first parameter)
#' # flat prior on baseline scale parameter (second parameter)
#' # and N(0,0.25) on regression parameter (third parameter)
#' 
#' inform_prior_fit <- ic_bayes(cbind(l, u) ~ grp, 
#'                              data = miceData,
#'                              logPriorFxn = priorFxn)
#' 
#' summary(flat_prior_model)
#' summary(inform_prior_fit)
#' # Note tight prior on the regression pulls posterior mean toward 0
#' 
#' @author Clifford Anderson-Bergman
#' @export
ic_bayes <- function(formula, data, logPriorFxn = function(x) return(0),
                     model = 'ph', dist = 'weibull', weights = NULL,
                     controls = bayesControls(), useMCores = F){

  if(missing(data)) data <- environment(formula)
  checkFor_cluster(formula)

  # Extracting x,y from formula + data
  reg_items = make_xy(formula, data)
  x = reg_items$x
  yMat = reg_items$y
  xNames = reg_items$xNames

  testMat <- cbind(x, 1)
  invertResult <- try(diag(solve(t(testMat) %*% testMat )), silent = TRUE)
  if(is(invertResult, 'try-error') & controls$useMLE_start){
    errorMsg <- 'covariate matrix is computationally singular!'
    errorMsg <- paste0(errorMsg, '\nic_bayes can still work with informative priors if useMLE_start = F in controls.')
    errorMsg <- paste0(errorMsg, '\nSee ?bayesControls for more details')
    stop( errorMsg )
  }
  
  modelName <- paste(dist, model, 'bayes')
  callText <- match.call()
  
  if(is.null(weights))	weights = rep(1, nrow(yMat))
  if(length(weights) != nrow(yMat))	stop('weights improper length!')
  if(min(weights) < 0)				stop('negative weights not allowed!')
  if(sum(is.na(weights)) > 0)			stop('cannot have weights = NA')

  # Recentering covariates
  covarOffset <- icColMeans(x)
  x <- t(t(x) - covarOffset)
  
  ans <- fit_bayes(yMat, x, parFam = dist, link = model, 
                   leftCen = 0, rightCen = Inf, uncenTol = 10^-6, 
                   regnames = xNames, weights = weights,
                   callText = callText, logPriorFxn = logPriorFxn,
                   bayesList = controls, modelName = modelName,
                   chains = controls$chains, use_m_cores = useMCores)
  ans$model = model
  
  # Extracting info needed for re-expanding data
  call_base = match.call(expand.dots = FALSE)
  call_info = readingCall(call_base)
  ans$terms <- call_info$mt
  ans$xlevels <- .getXlevels(call_info$mt, call_info$mf)
  ans$formula <- formula
  dataEnv <- new.env()
  dataEnv$data <- data
  ans$.dataEnv <- dataEnv
  ans$covarOffset <- matrix(covarOffset, nrow = 1)
  return(ans)
}

#' Control parameters for ic_bayes
#' 
#' @param samples               Number of samples. 
#' @param chains                Number of MCMC chains to run
#' @param useMLE_start          Should MLE used for starting point? 
#' @param burnIn                Number of samples discarded for burn in
#' @param samplesPerUpdate      Number of iterations between updates of proposal covariance matrix
#' @param initSD                If \code{useMLE_start == FALSE}, initial standard deviation used 
#' @param updateChol            Should cholesky decomposition be updated?
#' @param acceptRate            Target acceptance rate
#' @param thin                  Amount of thinning
#' 
#' @details 
#' 
#' Control parameters for the MH block updater used by \code{ic_bayes}.
#' 
#' The \code{samples} argument dictates how many MCMC samples are taken. One 
#' sample will be saved every \code{thin} iterations, so there will a total of
#' \code{thin * samples + burnIn} iterations. The burn in samples are not saved at all. 
#' 
#' Default behavior is to first calculate the MLE (not the MAP) estimate and use 
#' Hessian at the MLE to seed the proposal covariance matrix. After this, an updative 
#' covariance matrix is used. In cases with weakly informative likelihoods, 
#' using the MLE startpoint may lead to overly diffuse proposal or even undefined 
#' starting values. 
#' In this case, it suggested to use a cold start by setting \code{useMLE_start = F}
#' for the \code{controls} argument. In this case, the initial starting proposal
#'  covariance matrix will be a diagonal matrix with \code{initSD} standard deviations. 
#' 
#' @export
bayesControls <- function(samples = 4000, chains = 4,
                          useMLE_start = TRUE, burnIn = 2000, 
                          samplesPerUpdate = 1000, initSD = 0.1,
                          updateChol = TRUE, acceptRate = 0.25,
                          thin = 5){
  ans <- list(useMLE_start        = useMLE_start,
              chains              = chains,
              samples             = samples,
              thin                = thin,
              samplesPerUpdate    = samplesPerUpdate,
              updateChol          = updateChol,
              initSD              = initSD,
              burnIn              = burnIn,
              acceptRate          = acceptRate
  )
  return(ans)
}


fit_bayes <- function(y_mat, x_mat, parFam, link, 
                      leftCen = 0, rightCen = Inf, 
                      uncenTol = 10^-6, regnames, 
                      weights, callText, logPriorFxn,
                      bayesList, modelName, 
                      chains = 4, use_m_cores = T){
  parList<- make_par_fitList(y_mat, x_mat, parFam = parFam, 
                             link = link, leftCen = 0, rightCen = Inf,
                             uncenTol = 10^-6, regnames = regnames,
                             weights = weights, callText = modelName)
  bayesList$samples <- ceiling(bayesList$samples / chains)
  `%myDo%` <- `%do%`
  if(use_m_cores) `%myDo%` <- `%dopar%`
  seeds = runif(chains, 1, 10000)
  
  # fooling CRAN check...because it gets fooled by foreach
  this_seed = NULL
  
  c_fit_list <- foreach(this_seed = seeds) %myDo% {
    set.seed(this_seed)
    R_ic_bayes(bayesList, logPriorFxn, parList)
  }
  c_fit = c_fit_list[[1]]
  mcmc_samples = mcmc.list()
  allParNames = c(parList$bnames, parList$regnames)
  logPostDens <- list()
  for(i in 1:chains){
    these_samps <- c_fit_list[[i]]$samples
    colnames(these_samps) = allParNames
    these_samps <- coda::mcmc(these_samps, 
                              thin = bayesList$thin, 
                              start = bayesList$burnIn + 1)
    mcmc_samples[[i]] = these_samps
    logPostDens[[i]]  = c_fit_list[[i]]$logPosteriorDensity
  }

  nBase             <- length(parList$bnames)
  nRegPar           <- length(parList$regnames)
  
  fit <- new(modelName)
  fit$par           <- parFam
  mat_samples       <- as.matrix(mcmc_samples)
  fit$baseline      <- icr_colMeans(mat_samples[ ,1:nBase])
  regParVec         <- NULL
  covMat            <- NULL
  if(nRegPar > 0) { 
    regParVec <-  icr_colMeans(mat_samples[ ,nBase + 1:nRegPar])
  } 
  fit$reg_pars      <- regParVec
  fit$nSamples      <- nrow(mat_samples)
  fit$var           <- cov(mat_samples)
  fit$mcmcList      <- mcmc_samples
  fit$logPosteriorDensities <- logPostDens
  fit$ess           <- coda::effectiveSize(mcmc_samples)
  fit$call          <- callText
  fit$logPrior      <- logPriorFxn
  fit$finalChol     <- c_fit$finalChol
  MAP_info          <- getMaxPostDensInfo(logPostDens)
  MAP_dens          <- MAP_info[1]
  MAP_ind           <- MAP_info[2]
  fit$MAP_ind       <- MAP_ind
  fit$MAP_dens      <- MAP_dens
  fit$MAP_reg_pars  <- mat_samples[fit$MAP_ind, nBase + seq_len(nRegPar)]
  fit$MAP_baseline  <- mat_samples[fit$MAP_ind, 1:nBase]
  fit$samples   <- mat_samples
  if(nRegPar > 0) { 
    names(fit$reg_pars)    <- parList$regnames
    names(fit$MAP_reg_pars)    <- parList$regnames
  }
  names(fit$baseline)    <- parList$bnames
  names(fit$MAP_baseline)    <- parList$bnames
  fit$coefficients <- c(fit$baseline, fit$reg_pars)
  return(fit)
}


getMaxPostDensInfo <- function(maxDensList){
  vals <- NULL
  for(i in seq_along(maxDensList)){
    these_dens <- maxDensList[[i]]
    vals <- c(vals, these_dens)
  }  
  maxVal <- max(vals)
  maxValInd <- which(maxVal == maxVal)[1]
  ans <- c(maxVal, maxValInd)
  names(ans) <- c('maxVal', 'maxValInd')
  return(ans)
}

postDens_at_mle <- function(mle_fit, priorFxn){
  mle_vals <- coef(mle_fit)
  llk <- mle_fit$llk
  prior_dens <- priorFxn(mle_vals)
  post_dens <- llk + prior_dens
  ans <- c(post_dens, mle_vals)
  return(ans)
}