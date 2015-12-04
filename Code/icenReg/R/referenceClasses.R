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
                                  'final_llk',
                                  'iterations',
                                  'score'
                       ),
                       methods = list(
                         show = function(){
                           print(summary(.self) ) 
                         }
                       )
                       )

sp_fit <- setRefClass(Class = 'sp_fit',
                      contains = 'icenReg_fit',
                      fields = c('p_hat',
                                 'T_bull_Intervals',
                                 'bsMat', 
                                 'coef_bc'))

ic_ph <- setRefClass(Class = 'ic_ph',
                     contains = 'sp_fit')

ic_po <- setRefClass(Class = 'ic_po',
                     contains = 'sp_fit')

par_fit <- setRefClass(Class = 'par_fit',
                       contains = 'icenReg_fit',
                       fields = c('baseline',
                                  'hessian',
                                  'pca_hessian'
                                  ))

surv_trans_models <- c('po', 'ph')
parametricFamilies <- c('exponential', 'weibull', 'gamma', 'lnorm', 'loglogistic', 'generalgamma')

for(mod in surv_trans_models){
  for(fam in parametricFamilies)
    setRefClass(Class = paste(fam, mod, sep = " "),
                contains = 'par_fit')
}


###			Summary Class

setRefClass('icenRegSummary',
            fields = c('summaryParameters',
                       'model', 
                       'call', 
                       'baseline',
                       'sigFigs',
                       'fullFit',
                       'final_llk',
                       'iterations',
                       'other'),
            methods = list(
              initialize = function(fit){
                sigFigs <<- 4
                fullFit <<- fit
                model  <<- if(fit$model == 'ph') 'Cox PH' else 'Proportional Odds'
                baseline <<- fit$par
                colNames <- c('Estimate', 'Exp(Est)', 'Std.Error', 'z-value', 'p')
                coefs <- fit$coefficients
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
                summaryParameters <<- sumPars
                call <<- fit$call
                final_llk <<- fit$final_llk
                iterations <<- fit$iterations
                otherList <- list()
                if(inherits(fit, 'sp_fit')){
                  otherList[['bs_samps']] <- max(c(nrow(fit$bsMat),0))
                }
                other <<- otherList
              },
              show = function(){
                printSE <- TRUE
                sampSizeWarn <- FALSE
                if(baseline == 'semi-parametric'){
                  if(other[['bs_samps']] <= 1) printSE <- FALSE
                  if(other[['bs_samps']] < 100) sampSizeWarn <- TRUE
                }
                cat("\nModel: ", model, "\nBaseline: ", baseline, "\nCall: ")
                print(call)
                cat('\n')
                printMat <- summaryParameters
                if(!printSE) printMat <- printMat[,1:2]	
                print(printMat)
                cat('\nfinal llk = ', final_llk, '\nIterations = ', iterations, '\n')
                if(inherits(fullFit, 'sp_fit')) cat('Bootstrap Samples = ', other[['bs_samps']], '\n')
                if(sampSizeWarn){
                  cat("WARNING: only ", other[['bs_samps']], " bootstrap samples used for standard errors. Suggest using more bootstrap samples for inference\n")
                }
              }
            )
)
