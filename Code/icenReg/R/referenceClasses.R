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
                                  'pca_coefs',
                                  'pca_info',
                                  'baseOffset',
                                  'llk',
                                  'iterations',
                                  'score'
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
                         fields = c('samples',
                                    'baseline',
                                    'logPosteriorDensities',
                                    'nSamples',
                                    'ess',
                                    'logPrior',
                                    'finalChol', 
                                    'MAP_ind',
                                    'MAP_reg_pars',
                                    'MAP_baseline'))


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
                  sumPars <- summary(fit$samples)
                  otherList[['MAP']] <- signif( fullFit$samples[fullFit$MAP_ind,], sigFigs)
#                  names(otherList[['MAP']]) <- colnames(fullFit$samples)
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
                         )
                         )