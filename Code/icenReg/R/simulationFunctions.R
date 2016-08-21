#' Simulates Event Time from Survival Regression Model
#' 
#' @description Simulates event times from proportional hazards,
#' proportional odds and accelerated failure time models. 
#' 
#' @param linPred    A numeric vector of linear predictors for
#' simulated data
#' @param model      Regression model type. Options are \code{"ph"}, 
#' \code{"po"} or \code{"aft"} 
#' @param dist       The baseline distrubtion \code{q} function
#' @param paramList  A list of parameters to be passed to the baseline
#' distribution function 
simEventTime <- function(linPred = 0, model = 'ph', 
                         dist = qweibull, 
                         paramList = list(shape = 1, scale = 1)){
  n = length(linPred)
  rawP <- runif(n)
  nu <- exp(linPred)
  if(model == 'aft'){
    paramList$p <- rawP
    rawTimes <- do.call(dist, paramList)
    ans <- rawTimes * nu
    return(ans)
    
  }
  if (model == "ph") 
    adjFun <- function(x, nu) {
      1 - x^(1/nu)
    }
  else if (model == "po") 
    adjFun <- function(x, nu) {
      1 - x * (1/nu)/(x * 1/nu - x + 1)
    }
  adjP <- adjFun(rawP, nu)
  paramList$p <- adjP
  ans <- do.call(dist, paramList)
  return(ans)
}


simIC_weib <- function (n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, 
          scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL, 
          prob_cen = 1) 
{
  rawQ <- runif(n)
  x1 <- runif(n, -1, 1)
  x2 <- 1 - 2 * rbinom(n, 1, 0.5)
  nu <- exp(x1 * b1 + x2 * b2)
  if (model == "ph") 
    adjFun <- function(x, nu) {
      1 - x^(1/nu)
    }
  else if (model == "po") 
    adjFun <- function(x, nu) {
      1 - x * (1/nu)/(x * 1/nu - x + 1)
    }
  adjQ <- adjFun(rawQ, nu)
  trueTimes <- qweibull(adjQ, shape = shape, scale = scale)
  obsTimes <- runif(n = n, max = inspectLength)
  if (!is.null(rndDigits)) 
    obsTimes <- round(obsTimes, rndDigits)
  l <- rep(0, n)
  u <- rep(0, n)
  caught <- trueTimes < obsTimes
  u[caught] <- obsTimes[caught]
  l[!caught] <- obsTimes[!caught]
  if (inspections > 1) {
    for (i in 2:inspections) {
      oldObsTimes <- obsTimes
      obsTimes <- oldObsTimes + runif(n, max = inspectLength)
      if (!is.null(rndDigits)) 
        obsTimes <- round(obsTimes, rndDigits)
      caught <- trueTimes >= oldObsTimes & trueTimes < 
        obsTimes
      needsCatch <- trueTimes > obsTimes
      u[caught] <- obsTimes[caught]
      l[needsCatch] <- obsTimes[needsCatch]
    }
  }
  else {
    needsCatch <- !caught
  }
  u[needsCatch] <- Inf
  if (sum(l > u) > 0) 
    stop("warning: l > u! Bug in code")
  isCensored <- rbinom(n = n, size = 1, prob = prob_cen) == 
    1
  l[!isCensored] <- trueTimes[!isCensored]
  u[!isCensored] <- trueTimes[!isCensored]
  if (sum(l == Inf) > 0) {
    allTimes <- c(l, u)
    allFiniteTimes <- allTimes[allTimes < Inf]
    maxFiniteTime <- max(allFiniteTimes)
    l[l == Inf] <- maxFiniteTime
  }
  return(data.frame(l = l, u = u, x1 = x1, x2 = x2))
}