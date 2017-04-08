plot.icenReg_fit <- function(x, y, fun = 'surv', 
                             plot_legend = T,
                             cis = T, 
                             lgdLocation = 'topright', 
                             xlab = "time", ...){
  if(inherits(x, 'impute_par_icph'))	stop('plot currently not supported for imputation model')
  argList <- list(...)
  colors <- argList$col
  if(missing(y)) y <- argList$newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- icr_nrow(newdata)
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
    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
    
    addList <- list(xlab = xlab, ylab = yName, 
                    xlim = range(as.numeric(ranges), finite = TRUE), 
                    ylim = c(0,1))
    firstPlotList<- addListIfMissing(addList, firstPlotList)
    do.call(plot, firstPlotList)
  }
  
  lines(x, newdata, cis = cis, ...)
  if(nRows > 1 & plot_legend){
    grpNames <- rownames(newdata)
    legend(lgdLocation, legend = grpNames, lwd = rep(1, length(grpNames) ), col = colors)
  }
}

lines.surv_cis <- function(x, y,...){
  if(missing(y)) y <- NULL
  x$all_lines(cols = y, ...)
}

lines.icenReg_fit <- function(x, y, fun = 'surv', cis = F, ...){
  argList <- list(...)
  colors <- argList$col
  if(missing(y)) y <- argList$newdata	
  newdata <- y
  nRows <- 1
  if(!is.null(newdata)) nRows <- icr_nrow(newdata)
  if(nRows > 1){
    if(length(colors) == 0) colors = 1:nRows
    use_colors = colors
    if(length(use_colors) == 1) use_colors = rep(use_colors, nRows)
    if(length(use_colors) != nRows) stop("number of length(col) must be 0, 1 or equal to number of rows of new data")
    for(i in 1:nRows){
      this_row <- get_dataframe_row(y, i)
      this_col <- use_colors[i]
      argList$x <- x
      argList$y <- this_row
      argList$col <- this_col
      do.call(lines.icenReg_fit, argList)
    }
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
    ranges <- matrix(nrow = nRows, ncol = 2)
    ranges[,1] <- getFitEsts(x, newdata = newdata, p = 0.05 )
    ranges[,2] <- getFitEsts(x, newdata = newdata, p = 0.95 )
    if(is.null(colors))  colors <- 1:nRows
    for(i in 1:icr_nrow(ranges)){
      grid = ranges[i,1] + 0:100/100 * (ranges[i,2] - ranges[i,1])
      est.s <- 1 - getFitEsts(x, newdata = subsetData_ifNeeded(i, newdata), q = grid)
      argList[['x']] <- grid
      argList[['y']] <- s_trans(est.s)
      argList[['col']] <- colors[i]
      do.call(lines, argList)
    }
    if(cis){
      cis <- survCIs(x, newdata)
      lines(cis, colors)
    }
  }
}
