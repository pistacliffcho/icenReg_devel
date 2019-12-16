#' Plotting for icenReg Fits
#' @param x icenReg fit
#' @param y new data.frame
#' @param newdata new data.frame (ignored if \code{y} is included)
#' @param fun Function to be plotted. Options include \code{"surv"} or \code{"cdf"}
#' @param plot_legend Should legend be plotted?
#' @param cis Should confidence/credible interval be plotted?
#' @param ci_level Confidence/credible interval
#' @param survRange Range of survival curve to be plotted
#' @param evalPoints Number of evaluations of survival curve to be plotted. 
#' @param lgdLocation Location of legend; see \code{?legend} for options
#' @param xlab Label of x-axis
#' @param ... additional arguments to be passed to the base \code{plot} function
#' @details Plots survival function from either an \code{ic_np, ic_sp, ic_par} or \code{ic_bayes}
#' object. If \code{newdata} is \code{NULL}, the baseline distribution is plotted. Otherwise,
#'  \code{newdata} should be a \code{data.frame} with each row containing a set 
#'  covariates for which the fit will be plotted. If multiple rows are included, 
#'  the lines will be colored and a legend will be created using the rownames of \code{newdata}.
#'  
#'  For \code{ic_np} and \code{ic_sp}, the MLE is plotted with no intervals (at the time
#'  of writing this, there is no formula for standard errors of baseline distributions 
#'  for these methods). 
#' 
#'  For \code{ic_par} and \code{ic_bayes}, the output plotted is directly extracted from
#'  \code{survCIs}.
#'  
#'  If the argument \code{col} is provided, it will be used to color each 
#'  survival function in order of the rows provided in \code{newdata}.
#'  
#' @examples
#'  # Fitting mice data set
#'  data(miceData)
#'  miceFit <- ic_par(cbind(l, u) ~ grp, data = miceData) 
#'  
#'  # Creating covariates we want plotted
#'  newData <- data.frame(grp = c("ce", "ge"))
#'  # Naming rows for legend
#'  rownames(newData) <- c("Conventional", "Germ-Free")
#'  
#'  plot(miceFit, newdata = newData, 
#'       col = c('blue', 'orange'))
#' @export
plot.icenReg_fit <- function(x, y, newdata = NULL, fun = 'surv', 
                             plot_legend = T,
                             cis = T, ci_level = 0.9,
                             survRange = c(0.025, 1), 
                             evalPoints = 200,
                             lgdLocation = lgd_default(fun), 
                             xlab = "time", ...){
  if(inherits(x, 'impute_par_icph'))	stop('plot currently not supported for imputation model')
  argList <- list(...)
  if(missing(y)) y <- newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- icr_nrow(newdata)
  colors <- argList$col
  if(is.null(colors)) colors = 1:nRows
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){ s_trans <- function(x){1-x}; yName = 'F(t)' }
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')
  
  addList <- list(xlab = xlab, ylab = yName)
  argList <- addListIfMissing(addList, argList)
  firstPlotList <- argList
  firstPlotList[['type']] <- 'n'
  firstPlotList[['x']] <- 1
  firstPlotList[['y']] <- 1
  
  
  if(x$par == 'semi-parametric' | x$par == 'non-parametric'){
    curveInfo <- getSCurves(x, y)
    allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
    dummyx <- range(allx, finite = TRUE)
    dummyy <- c(0,1)
    firstPlotList[['xlim']] = dummyx
    firstPlotList[['ylim']] = dummyy
    do.call(plot, firstPlotList)
    
  }
  else if(inherits(x, 'par_fit') | inherits(x, 'bayes_fit')){
    pRange = 1 - survRange
    ranges <- matrix(nrow = nRows, ncol = 2)
    surv_est <- survCIs(x, newdata, p = pRange)
    xmin = Inf
    xmax = -Inf
    for(i in seq_along(surv_est)){
      this_est <- surv_est$cis[[i]][,3]
      xmin = min(c(xmin, this_est))
      xmax = max(c(xmax, this_est))
    }
    
    addList <- list(xlab = xlab, ylab = yName, 
                    xlim = c(xmin, xmax), 
                    ylim = c(0,1), col = NA)
    firstPlotList<- addListIfMissing(addList, firstPlotList)
    do.call(plot, firstPlotList)
  }
  
  newLinesList = list(...)
  newLinesList$x = x
  newLinesList$y = newdata
  newLinesList$cis = cis
  newLinesList$ci_level = ci_level
  newLinesList$col = colors
  newLinesList$evalPoints = evalPoints
  newLinesList$survRange = survRange
  newLinesList$fun = fun
  do.call(lines.icenReg_fit, newLinesList)
  if(nRows > 1 & plot_legend){
    grpNames <- rownames(newdata)
    legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = colors)
  }
}

lines.surv_cis <- function(x, col, include_cis = T, fun = "surv", ...){ 
  x$all_lines(col = col, include_cis = include_cis, fun, ...) }

#' Function for default legend location
#' @param fun_type Type of function to plot. Options are 'surv' or 'cdf'
#' @description Returns default locations for plot legends
#' @noRd
lgd_default = function(fun_type = "surv"){
  if(fun_type == "surv") ans = "topright"
  else if(fun_type == "cdf") ans = "topleft"
  else ans = NA
  return(ans)
}

#' @inherit plot.icenReg_fit
#' @export
lines.icenReg_fit <- function(x, y, newdata = NULL, 
                              fun = 'surv', cis = F, 
                              ci_level = 0.9, 
                              survRange = c(0.025, 1), 
                              evalPoints = 20, ...){
  # Renaming for clarity
  model <- x
  # Resolving newdata
  if(missing(y)) y <- newdata	
  newdata <- y
  
  argList <- list(...)
  # Resolving colors 
  colors <- argList$col
  nRows <- 1
  if(!is.null(newdata)) nRows <- icr_nrow(newdata)
  if(nRows > 1){
    if(length(colors) == 0) colors = 1:nRows
    use_colors = colors
    if(length(use_colors) == 1) use_colors = rep(use_colors, nRows)
    if(length(use_colors) != nRows) 
      stop("number of length(col) must be 0, 1 or equal to number of rows of new data")
  }
  else{
    if(length(colors) == 0) colors = 1
  }
  
  # Setting up survival or cdf 
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){ s_trans <- function(x){1-x}; yName = 'F(t)' }
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')
  
  # Semi/non-parametric models
  if(model$par == 'semi-parametric' | model$par == 'non-parametric'){
    argList <- addIfMissing('s', 'type', argList)
    curveInfo <- getSCurves(model, newdata)
    allx <- c(curveInfo$Tbull_ints[,1], curveInfo$Tbull_ints[,2])
    dummyx <- range(allx, finite = TRUE)
    dummyy <- c(0,1)
    x_l <- curveInfo$Tbull_ints[,1]
    x_u <- curveInfo$Tbull_ints[,2]
    k <- length(x_l)
    ss <- curveInfo$S_curves
    if(is.null(colors))  colors <- 1:length(ss)
    if(length(colors) == 1) colors <- rep(colors, length(ss)) 
    for(i in 1:length(ss)){
      argList[['x']] <- x_l
      argList[['y']] <- s_trans(ss[[i]])
      argList[['col']] <- colors[i]
      do.call(lines, argList)     
      argList[['x']] <- x_u
      do.call(lines, argList)     
      argList[['x']] <- c(x_l[k], x_u[k])
      argList[['y']] <- s_trans(c(ss[[i]][k], ss[[i]][k]))
      do.call(lines, argList)
    }
  }
  # Parametric/Bayes models
  else if(inherits(x, 'par_fit') | inherits(x, 'bayes_fit')){
    pRange = 1 - survRange
    p_eval = pRange[1] + (pRange[2] - pRange[1]) * (0:evalPoints) / evalPoints
    cis_est <- survCIs(x, newdata, p = p_eval)
    argList <- c(list(fun = fun), list(...))
    argList$col = NULL
    argList$x = cis_est
    argList$col = colors
    argList$include_cis = cis
    argList$evalPoints = NULL
    argList$survRange = NULL
    do.call(lines.surv_cis, argList)
  }
}


plot.sp_curves <- function(x, sname = 'baseline', xRange = NULL, ...){
  if(is.null(xRange))
    xRange <- range(c(x[[1]][,1], x[[1]][,2]), finite = TRUE)
  dotList <- list(...)
  addList <- list(xlim = xRange, ylim = c(0,1), x = NA)
  dotList <- addListIfMissing(addList, dotList)
  do.call(plot, dotList)
  lines(x, sname = sname, ...)
}

lines.ic_npList <- function(x, fitNames = NULL, fun = "surv", ...){
  if(is.null(fitNames)){
    fitNames <- names(x$scurves)
    lines(x, fitNames, fun = fun, ...)
    return()
  }
  dotList <- list(...)
  dotList$fun = fun
  cols <- dotList$col
  
  for(i in seq_along(fitNames)){
    thisName <- fitNames[i]
    dotList$col <- cols[i]
    dotList$x <- x$scurves[[thisName]]
    do.call(lines, dotList)
  }
}

plot.ic_npList <- function(x, fitNames = NULL, 
                           fun = "surv",
                           lgdLocation = lgd_default(fun),
                           plot_legend = T,
                           ... ){
  # Renaming for clarity
  npList <- x
  ### Setting up plot frame  ###
  if(fun == "surv"){ yLab = "S(t)" }
  else if(fun == "cdf"){ yLab = "F(t)" }
  addList <- list(xlim = npList$xRange,
                  ylim = c(0,1),
                  xlab = 'time', 
                  ylab = yLab,
                  x = NA)
  dotList <- list(...) 
  possibleNames <- names(npList$scurves)
  if(any(names(dotList) == 'newdata')){
    nd <- dotList$newdata
    if(is.data.frame(nd)){ fitNames <- as.character(nd[[1]])}
    else fitNames <- as.character(nd)
    dotList$newdata = NULL
  }
  
  if(!all(fitNames %in% possibleNames)){
    this_warn <- paste('Invalid group names. \nValid names:', 
                       possibleNames, 
                       '\nProvided names:', 
                       fitNames, collapse = " ")
    stop(this_warn)
  }
  dotList <- addListIfMissing(addList, dotList)
  do.call(plot, dotList)  
  
  ### Plotting lines ###
  grpNames <- names(npList$fitList)
  cols <- dotList$col
  if(is.null(cols)) cols = 2:(length(grpNames) + 1)
  if(length(cols) != length(grpNames)) 
    stop('length of number of strata not equal to number of colors')
  dotList$col <- cols
  dotList$fitNames = fitNames
  dotList$x <- npList
  dotList$fun = fun
  do.call(lines, dotList)
  if(plot_legend)
   legend(lgdLocation, legend = grpNames, col = cols, lty = 1)
}

lines.sp_curves <- function(x, sname = 'baseline', fun = "surv", ...){
  # Renaming for clarity
  sp_curves <- x
  
  # Drawing first lines
  firstTimeObs <- sp_curves[[1]][1,1]
  firstTimeAssume <- firstTimeObs
  if(fun == "surv"){ first_y = 1 }
  else if(fun == "cdf"){ first_y = 0 }
  else stop("fun type no recognized. Options = 'surv' or 'cdf'")
  if(firstTimeObs > 0)
    firstTimeAssume <- 0
  lines(c(firstTimeAssume, firstTimeObs), rep(first_y, 2), ...)
  
  # Drawing step function
  if(fun == "surv") y = sp_curves[[2]][[sname]]
  else if(fun == "cdf") y = 1 - sp_curves[[2]][[sname]]
  lines(sp_curves[[1]][,1], y, ..., type = 's')
  lines(sp_curves[[1]][,2], y, ..., type = 's')
  
  # Drawing last lines
  lastObs <- nrow(x[[1]])
  lastTimes <- sp_curves[[1]][lastObs,]
  if(lastTimes[2] == Inf) lastTimes[2] <- lastTimes[1]
  lastTimes[2] <- lastTimes[2] + (lastTimes[2] - firstTimeObs)
  lines(lastTimes, 1 - rep(first_y, 2), ... ) 
}
