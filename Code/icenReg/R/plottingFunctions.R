plot.icenReg_fit <- function(x, y, newdata = NULL, fun = 'surv', 
                             plot_legend = T,
                             cis = T, ci_level = 0.9,
                             lgdLocation = 'topright', 
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
    ranges <- matrix(nrow = nRows, ncol = 2)
    surv_est <- survCIs(x, newdata, p = c(0.025, 0.975))
    xmin = Inf
    xmax = -Inf
    for(i in seq_along(surv_est)){
      this_est <- surv_est$cis[[i]][,3]
      xmin = min(c(xmin, this_est))
      xmax = max(c(xmax, this_est))
    }
    # ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    # ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
    
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
  do.call(lines, newLinesList)
  # lines(x, newdata, cis = cis, ci_level = ci_level,...)
  if(nRows > 1 & plot_legend){
    grpNames <- rownames(newdata)
    legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = colors)
  }
}

lines.surv_cis <- function(x, col,...){ x$all_lines(col = col, ...) }

lines.icenReg_fit <- function(x, y, newdata = NULL, 
                              fun = 'surv', cis = F, ci_level = 0.9, ...){
  argList <- list(...)
  colors <- argList$col
  if(missing(y)) y <- newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- icr_nrow(newdata)
  if(nRows > 1){
    if(length(colors) == 0) colors = 1:nRows
    use_colors = colors
    if(length(use_colors) == 1) use_colors = rep(use_colors, nRows)
    if(length(use_colors) != nRows) stop("number of length(col) must be 0, 1 or equal to number of rows of new data")
  }
  if(fun == 'surv'){ s_trans <- function(x){x}; yName = 'S(t)'}
  else if(fun == 'cdf'){ s_trans <- function(x){1-x}; yName = 'F(t)' }
  else stop('"fun" option not recognized. Choices are "surv" or "cdf"')
  if(x$par == 'semi-parametric' | x$par == 'non-parametric'){
    argList <- addIfMissing('s', 'type', argList)
    curveInfo <- getSCurves(x, y)
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
  else if(inherits(x, 'par_fit') | inherits(x, 'bayes_fit')){
    cis_est <- survCIs(x, newdata)
    argList <- list(...)
    argList$col = NULL
    argList$x = cis_est
    argList$col = colors
    argList$include_cis = cis
    do.call(lines, argList)
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

lines.ic_npList <- function(x, fitNames = NULL, ...){
  if(is.null(fitNames)){
    fitNames <- names(x$scurves)
    lines(x, fitNames, ...)
  }
  dotList <- list(...)
  cols <- dotList$col
  
  for(i in seq_along(fitNames)){
    thisName <- fitNames[i]
    dotList$col <- cols[i]
    dotList$x <- x$scurves[[thisName]]
    do.call(lines, dotList)
  }
}

plot.ic_npList <- function(x, fitNames = NULL, lgdLocation = 'bottomleft', ... ){
  addList <- list(xlim = x$xRange,
                  ylim = c(0,1),
                  xlab = 'time', 
                  ylab = 'S(t)', 
                  x = NA)
  dotList <- list(...)
  possibleNames <- names(x$scurves)
  if(any(names(dotList) == 'newdata')){
    nd <- dotList$newdata
    if(is.data.frame(nd)){ fitNames <- as.character(nd[[1]])}
    else fitNames <- as.character(nd)
    dotList$newdata = NULL
  }
  
  if(!all(fitNames %in% possibleNames)){
    this_warn <- paste('Invalid group names. \nValid names:', possibleNames, 
                       '\nProvided names:', fitNames, collapse = " ")
    stop(this_warn)
  }
  dotList <- addListIfMissing(addList, dotList)
  do.call(plot, dotList)  
  grpNames <- names(x$fitList)
  cols <- dotList$col
  if(is.null(cols)) cols = 2:(length(grpNames) + 1)
  if(length(cols) != length(grpNames)) 
    stop('length of number of strata not equal to number of colors')
  dotList$col <- cols
  dotList$fitNames = fitNames
  dotList$x <- x
  do.call(lines, dotList)
  legend(lgdLocation, legend = grpNames, col = cols, lty = 1)
}

lines.sp_curves <- function(x, sname = 'baseline',...){
  firstTimeObs <- x[[1]][1,1]
  firstTimeAssume <- firstTimeObs
  if(firstTimeObs > 0)
    firstTimeAssume <- 0
  lines(c(firstTimeAssume, firstTimeObs), c(1,1), ...)
  lines(x[[1]][,1], x[[2]][[sname]], ..., type = 's')
  lines(x[[1]][,2], x[[2]][[sname]],   ..., type = 's')
  lastObs <- nrow(x[[1]])
  lastTimes <- x[[1]][lastObs,]
  if(lastTimes[2] == Inf) lastTimes[2] <- lastTimes[1]
  lastTimes[2] <- lastTimes[2] + (lastTimes[2] - firstTimeObs)
  lines(lastTimes, c(0,0), ... ) 
}
