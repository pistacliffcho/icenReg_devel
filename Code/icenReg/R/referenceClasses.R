icenReg_fit <- setRefClass(Class = 'icenReg_fit',
                       fields = c('var',
                                  'coefficients',
                                  'formula',
                                  'call',
                                  '.dataEnv',
                                  'par',
                                  'model',
                                  'reg_pars',
                                  'terms',
                                  'xlevels',
                                  'llk',
                                  'iterations',
                                  'score', 
                                  'covarOffset', 
                                  'depType'
                       ),
                       methods = list(
                         show = function(){
                           print(summary(.self) ) 
                         },
                         getRawData = function(){
                           return(.self$.dataEnv$data)
                         },
                         initialize = function(...){
                           makeActiveBinding('coef', function(){return(coefficients)} ,.self)
                           depType <<- "Independence"
                         }
                       )
                       )

sp_fit_class <- setRefClass(Class = 'sp_fit',
                      contains = 'icenReg_fit',
                      fields = c('p_hat',
                                 'T_bull_Intervals',
                                 'bsMat', 
                                 'coef_bc'))

ic_np_class <- setRefClass(Class = 'ic_np', 
                     contains = 'sp_fit')

ic_ph_class <- setRefClass(Class = 'ic_ph',
                     contains = 'sp_fit')

ic_po_class <- setRefClass(Class = 'ic_po',
                     contains = 'sp_fit')


par_class <- setRefClass(Class = 'par_fit',
                       contains = 'icenReg_fit',
                       fields = c('baseline',
                                  'hessian',
                                  'pca_hessian'
                                  ))

bayes_fit <- setRefClass(Class = 'bayes_fit',
                         contains = 'icenReg_fit',
                         fields = c('mcmcList',
                                    'baseline',
                                    'logPosteriorDensities',
                                    'nSamples',
                                    'ess',
                                    'logPrior',
                                    'finalChol', 
                                    'MAP_dens',
                                    'MAP_ind',
                                    'MAP_reg_pars',
                                    'MAP_baseline', 
                                    'samples'))


surv_trans_models <- c('po', 'ph', 'aft', 'none')
parametricFamilies <- c('exponential', 'weibull', 'gamma', 'lnorm', 'loglogistic', 'generalgamma')

for(mod in surv_trans_models){
  for(fam in parametricFamilies){
    setRefClass(Class = paste(fam, mod, sep = " "),
                contains = 'par_fit')
    setRefClass(Class = paste(fam, mod, 'bayes', sep = " "),
              contains = 'bayes_fit')
  }  
}


###			Summary Class

setRefClass('icenRegSummary',
            fields = c('summaryParameters',
                       'model', 
                       'call', 
                       'baseline',
                       'sigFigs',
                       'fullFit',
                       'llk',
                       'iterations',
                       'other'),
            methods = list(
              initialize = function(fit){
                sigFigs <<- 4
                fullFit <<- fit
                otherList <- list()
                if(fit$model == 'ph') model <<- 'Cox PH'
                if(fit$model == 'po') model <<- 'Proportional Odds'
                if(fit$model == 'aft') model <<- 'Accelerated Failure Time'
                if(fit$model == 'none') model <<- 'Non-parametric'
                if(inherits(fit, 'bayes_fit')) model <<- paste("Bayesian", model)
                baseline <<- fit$par
                colNames <- c('Estimate', 'Exp(Est)', 'Std.Error', 'z-value', 'p')
                coefs <- fit$coefficients
                if(is(fit, 'bayes_fit')){
                  sumPars <- summary(fit$mcmcList)
                  otherList[['MAP']] <- signif( fullFit$samples[fullFit$MAP_ind,], sigFigs)
                  # otherList[["Gelman"]] <- coda::gelman.diag(fit$mcmcList)
                }
                else{
                  sumPars <- matrix(nrow = length(coefs), ncol = length(colNames))
                  se <- sqrt(diag(fit$var))
                  for(i in seq_along(coefs)){
                    sumPars[i,1] <- coefs[i]
                    sumPars[i,2] <- exp(coefs[i])
                    sumPars[i,3] <- se[i]
                    sumPars[i,4] <- coefs[i]/se[i]
                    sumPars[i,5] <- 2 * (1 - pnorm(abs(sumPars[i,4])))
                  }
                  colnames(sumPars) <- colNames
                  rownames(sumPars) <- names(coefs)
                  sumPars <- signif(sumPars, sigFigs)
                }
                summaryParameters <<- sumPars
                call <<- fit$call
                llk <<- fit$llk
                iterations <<- fit$iterations
                if(inherits(fit, 'sp_fit') & !inherits(fit, 'ic_np')){
                  otherList[['bs_samps']] <- max(c(nrow(fit$bsMat),0))
                }
                other <<- otherList
              },
              show = function(){
                printSE <- TRUE
                sampSizeWarn <- FALSE
                if(inherits(fullFit, 'ic_np')) printSE = FALSE
                if(baseline == 'semi-parametric'){
                  if(other[['bs_samps']] <= 1) printSE <- FALSE
                  if(other[['bs_samps']] < 100) sampSizeWarn <- TRUE
                }
                cat("\nModel: ", model)
                cat("\nDependency structure assumed:", fullFit$depType)
                if(!inherits(fullFit, 'ic_np')){
                  cat("\nBaseline: ", baseline, "\nCall: ")
                  print(call)
                }
                cat('\n')
                printMat <- summaryParameters
                if(!printSE){
                  if(nrow(printMat) > 1) printMat <- printMat[,1:2]
                  else{
                    if(nrow(printMat) == 0) {printMat <- "No covariates used \n"}
                    else{
                      newMat <- matrix(nrow = 1, ncol = 2)
                      newMat[1,] <- printMat[1,1:2]
                      rownames(newMat) <- rownames(printMat)
                      colnames(newMat) <- colnames(printMat)[1:2]
                      printMat <- newMat
                    }
                  }
                }
                if(is.character(printMat)) { cat(printMat)}
                else{print(printMat)}
                if(!inherits(fullFit, 'bayes_fit')) cat('\nfinal llk = ', llk, '\nIterations = ', iterations, '\n')
                if(inherits(fullFit, 'sp_fit') & !inherits(fullFit, 'ic_np')) cat('Bootstrap Samples = ', other[['bs_samps']], '\n')
                if(inherits(fullFit, 'bayes_fit')){
                  cat('3. MAP estimates:\n') 
                  print( other[['MAP']] )
                  # cat('\n4. Gelman Diagnostics:\n')
                  # print(other[['Gelman']])
                }
                if(sampSizeWarn){
                  cat("WARNING: only ", other[['bs_samps']], " bootstrap samples used for standard errors. \nSuggest using more bootstrap samples for inference\n")
                }
              }
            )
)


ic_npList <- setRefClass(Class = 'ic_npList',
                         fields = c('fitList', 'xRange', 'scurves', 'nGrp'),
                         methods = list(
                           initialize = function(fitList){
                             fitList <<- fitList
                             xVals <- c(Inf, -Inf)
                             scList <- list()
                             grpCounts <- numeric()
                             for(fitName in names(fitList)){
                               fit <- fitList[[fitName]]
                               thisSC <- getSCurves(fit)
                               scList[[fitName]] <- thisSC
                               xVals <- range(c(thisSC$Tbull_ints, xVals), finite = TRUE )
                            #   xVals[1] <- min( c(thisSC$Tbull_ints[1], xVals[1]) )
                            #   xVals[2] <- max( c(tail(thisSC$Tbull_ints[,2], 1), xVals[2]) )
                               grpCounts[fitName] <- nrow(getData(fit))
                             }
                             nGrp <<- grpCounts
                             xRange <<- xVals
                             scurves <<- scList
                           },
                           show = function(){
                             cat("Stratified NPMLE for interval censored data")
                             cat("\nGroup Counts:\n")
                             print(nGrp)
                           }
                         ))


surv_cis <- setRefClass("surv_cis",
                          fields = c('cis', 'call', 'newdata', 'ci_level'),
                          methods = list(
                            fit_one_row = function(fit, newdata_row,
                                                   p, q, p_ends, MC_samps){
                              samp <- sampleSurv(fit, newdata_row, 
                                               p = p, q = q, 
                                               samples = MC_samps)
                              
                              use_input <- 'p'
                              inputLength <- length(p)
                              inputVals <- p
                              if(is.null(p)){
                                use_input = 'q'
                                inputLength = length(q)
                                inputVals <- q
                              }
                              q_ests <- matrix(nrow = inputLength, ncol = 3)
                              mean_ests <- NULL 
                              p_use <- c(0.5, p_ends)
                              for(i in 1:inputLength){ 
                                q_ests[i,] <- quantile(samp[,i], probs = p_use) 
                                mean_ests[i] <- mean(samp[,i])
                              }
                              ans <- cbind(inputVals, mean_ests, q_ests)
                              firstColName <- 'Percentile'
                              if(is.null(p)) firstColName <- 'Time'
                              colnames(ans) <- c(firstColName, 
                                                 'estimate (mean)', 'estimate (median)', 
                                                 "lower", "upper")
                              return(ans)
                            },
                            initialize = function(fit, 
                                                  newdata = NULL,
                                                  p = c(0:9 * .1 + 0.05),
                                                  q = NULL, 
                                                  ci_level = 0.95, 
                                                  MC_samps = 40000){
                              call <<- fit$call
                              ci_level <<- ci_level
                              if(ci_level < 0 | ci_level > 1) stop('invalid ci_level')
                              alpha <- (1 - ci_level)/2
                              p_low = alpha
                              p_hi  = 1 - alpha
                              ci_list <- list()
                              rowNames <- rownames(newdata)
                              if(is.null(newdata)){
                                rowNames <- 'baseline'
                              }
                              else if(is.null(rowNames)){
                                if(nrow(newdata) > 0){
                                  rownames(newdata) <<- 1:nrow(newdata)
                                  rowNames <- rownames(newdata)
                                }
                              }
                              for(i in seq_along(rowNames)){
                                this_name <- rowNames[i]
                                ci_list[[this_name]] <- fit_one_row(fit, newdata_row = get_dataframe_row(newdata, i), 
                                                                    p = p, q = q, 
                                                                    p_ends = c(p_low, p_hi), 
                                                                    MC_samps = MC_samps)
                              }
                              cis <<- ci_list
                            },
                            show = function(){
                              cat("Model call:\n  ")
                              print(call,) 
                              cat("Credible Level =", ci_level, "\n")
                              for(i in seq_along(cis)){
                                this_name <- names(cis)[i]
                                cat("Rowname: ", this_name, "\n")
                                print(cis[[i]])
                              }
                            },
                            one_lines = function(index = 1, this_col, 
                                                 include_cis, fun = "surv", ...){
                              argList <- list(...)
                              argList$col = this_col
                              argList$lty = 1
                              these_cis <- cis[[index]]
                              perc = these_cis[,1]
                              est = these_cis[,3]
                              lower = these_cis[,4]
                              upper = these_cis[,5]
                              if(fun == "surv"){
                                argList$y = 1 - perc
                              }
                              else if(fun == "cdf"){
                                argList$y = perc
                              }
                              argList$x = est
                              do.call(lines, argList)
                              if(include_cis){
                                argList$lty = 2
                                argList$x = lower
                                do.call(lines, argList)
                                argList$x = upper
                                do.call(lines, argList)
                              }
                            },
                            all_lines = function(cols = NULL, include_cis,
                                                 fun = "surv", ...){
                              nCIs = length(cis)
                              nCols <- length(cols)
                              if(nCols == 0){ cols = 1:nCIs; nCols = length(cols) }
                              if(nCols == 1){ cols = rep(cols, nCIs); nCols = length(cols) }
                              if(nCols != nCIs) stop("number colors provided does not match up with number of CI's to plot")
                              argList <- c(list(fun = fun), list(...))
                              argList$include_cis = include_cis
                              for(i in seq_along(cis)){ 
                                argList$index = i
                                argList$this_col = cols[i]
                                argList$cols <- NULL
                                do.call(.self$one_lines, argList)
                                }
                            }
                          ))