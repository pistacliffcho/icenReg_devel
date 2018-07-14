#' Non-Parametric Estimator for Interval Censored Data
#'
#' @description
#'  Fits the non-parametric maximum likelihood estimator (NPMLE) for univariate interval censored data. 
#'  This is a generalization of the Kaplan-Meier curves that allows for interval censoring. 
#'  Also referred to as the Turnbull estimator.
#' 
#' @param formula   Formula for stratification. If only one group, can be left blank and 
#' data must be entered as n x 2 matrix.
#' @param data      A \code{data.frame} or an n x 2 matrix. See details.
#' @param maxIter   Maximum iterations
#' @param tol       Numeric tolerance
#' @param B         Should intervals be open or closed? See details.
#' @param weights   Weights (optional)
#'  
#' @details 
#' \code{data} must be an n x 2 matrix or data.frame containing two columns of data representing 
#' left and right sides of the censoring interval, denoted L and R. This allows for left censored 
#' (L == 0), right censored (R == inf), uncensored (L == R) along with general interval censored observations. 
#'
#' The argument \code{B} determines whether the intervals should be open or closed, i.e. 
#' \code{B = c(0,1)} implies that the event occurs within the interval \code{(l,u]}. 
#' The exception is that if \code{l == u}, it is assumed that the event is uncensored, regardless of \code{B}.
#'
#' The NPMLE is fit using an efficient implementation of the EMICM algorithm.
#' @examples 
#' data(miceData)
#' fit <- ic_np(cbind(l, u) ~ grp, data = miceData)
#' # Stratifies fits by group
#' 
#' plot(fit) 
#' @references 
#' Turnbull, B. (1976) The empricial distribution with arbitrarily grouped and censored data 
#' \emph{Journal of the Royal Statistical Society B}, vol 38 p290-295
#'
#' Wellner, J. A., and Zhan, Y. (1997) A hybrid algorithm for computation of the maximum likelihood estimator 
#' from censored data, \emph{Journal of the  American Statistical Association}, Vol 92, pp945-959
#'
#' Anderson-Bergman, C. (2016) An efficient implementation of the EMICM algorithm for the interval censored NPMLE
#' \emph{Journal of Computational and Graphical Statistics}, \emph{just accepted}
#' 
#' @author Clifford Anderson-Bergman
#' @export
ic_np <- function(formula = NULL, data, maxIter = 1000, tol = 10^-10, B = c(0,1), 
                  weights = NULL){
  if(!inherits(formula, 'formula')) {
    #Covering when user ONLY provides data as first unlabeled argument
    data <- formula
    formula <- NULL
  }
  if(is.null(weights)) weights = rep(1, nrow(data))
  if(is.null(formula)){ 
    ans = ic_npSINGLE(data, maxIter = maxIter, tol = tol, B = B, weights = weights)
    return(ans) 
  }
  
  if(missing(data)) data <- environment(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  callInfo <- readingCall(mf)
  mf <- callInfo$mf
  mt <- callInfo$mt

  y <- model.response(mf, "numeric")
  yMat <- makeIntervals(y, mf)
  yMat <- adjustIntervals(B, yMat)
  
  
  formFactor <- formula[[3]]
  if( formula_has_plus(formFactor) ){ 
    stop('predictor must be either single factor OR 0 for ic_np')
  }
  if(formFactor == 0){ return(ic_npSINGLE(yMat, maxIter = maxIter, tol = tol, B = B, weights = weights)) }
  #  thisFactor <- data[[ as.character(formFactor) ]]
  thisFactor <- getFactorFromData(formFactor, data)
  if(!is.factor(thisFactor)){ stop('predictor must be factor') }
  theseLevels <- levels(thisFactor)
  fitList <- list()
  for(thisLevel in theseLevels){
    thisData <- yMat[thisFactor == thisLevel, ]
    this_w = weights[thisFactor == thisLevel]
    if(nrow(thisData) > 0)
      fitList[[thisLevel]] <- ic_npSINGLE(thisData, maxIter = maxIter, 
                                          tol = tol, B = B, 
                                          weights = this_w)
  }
  ans <- ic_npList(fitList)
  return(ans)
}

ic_npSINGLE <- function(data,  maxIter = 1000, tol = 10^-10, B, weights){
  data <- as.matrix(data)
  # Weights must be non-negative
  if(any(weights < 0)) stop("weights supplied cannot be less than 0")
  # Our algorithm requires positive weights
  zeroWeight = weights == 0
  data = data[!zeroWeight, ]
  weights = weights[!zeroWeight]
  
  if(ncol(data) != 2) stop("data should be an nx2 matrix or data.frame")
  if(any(data[,1] > data[,2]) ) stop(paste0("data[,1] > data[,2].",
                                            "This is impossible for interval censored data") )
  storage.mode(data) <- "double"
  data <- adjustIntervals(B, data)
  mis <- findMaximalIntersections(data[,1], data[,2])
  fit <- .Call("EMICM", mis$l_inds, mis$r_inds, as.integer(maxIter), weights)
  tbulls <- rbind(mis$mi_l, mis$mi_r)
  ans <- new('ic_np')
  #ans <- list(phat = fit[[1]], Tbull_ints = tbulls, llk = fit[[2]], iters = fit[[3]])
  ans$p_hat <- fit[[1]]
  ans$T_bull_Intervals <- tbulls
  ans$coefficients <- numeric()
  ans$llk <- fit[[2]]
  ans$iterations <- fit[[3]]
  ans$par <- 'non-parametric'
  ans$model = 'none'
  ans$var <- matrix(nrow = 0, ncol = 0) 
  dataEnv <- new.env()
  dataEnv[['data']] <- data
  dataEnv[["weights"]] = weights
  ans[['.dataEnv']] <- dataEnv
  return(ans)
}
