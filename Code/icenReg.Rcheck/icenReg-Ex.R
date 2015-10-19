pkgname <- "icenReg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "icenReg-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('icenReg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ICNPMLE")
### * ICNPMLE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ICNPMLE
### Title: Computes the NPMLE for Univariate or Bivariate Interval Censored
###   Data
### Aliases: ICNPMLE IC_NPMLE

### ** Examples

simData <- simBVCen(500)

fit <- ICNPMLE(simData)

fit



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ICNPMLE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diag_baseline")
### * diag_baseline

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diag_baseline
### Title: Compare parametric baseline distributions with semi-parametric
###   baseline
### Aliases: diag_baseline

### ** Examples

	 data(essIncData_small)
	 useData <- essIncData_small
	 
	 #Note to user: suggest replacing useData with essIncData
	 #instead of essIncData_small. Using small dataset to quickly 
	 #pass CRAN tests
	 
	 fit_po <- ic_sp(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry,
	                 data = useData, bs_samples = 0, model = 'po')
	
	 diag_baseline(fit_po)

	 #Weibull model appears best fit



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diag_baseline", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diag_covar")
### * diag_covar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diag_covar
### Title: Evaluate covariate effect for regression model
### Aliases: diag_covar

### ** Examples

	 data(essIncData_small)
	 useData <- essIncData_small
	 
	 #Note to user: suggest replacing useData with essIncData
	 #instead of essIncData_small. Using small dataset to quickly 
	 #pass CRAN tests
	 par(mfrow = c(1,2))
	 	 
	 diag_covar(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 
		         data = useData, model = 'po')
		
	 diag_covar(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 
		         data = useData, model = 'ph')


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diag_covar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("essIncData")
### * essIncData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: essIncData
### Title: Interval Censored Income Data from European Social Survey
### Aliases: essIncData

### ** Examples

	data(essIncData)
	
	lnormFit <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 	
	                   data = essIncData,
	                   model = 'po',
	                   dist = 'loglogistic')
	
	summary(lnormFit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("essIncData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("essIncData_small")
### * essIncData_small

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: essIncData_small
### Title: Interval Censored Income Data from European Social Survey
### Aliases: essIncData_small

### ** Examples

	data(essIncData_small)
	
	lnormFit <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 	
	                   data = essIncData_small,
	                   model = 'po',
	                   dist = 'loglogistic')
	
	summary(lnormFit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("essIncData_small", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("getFitEsts")
### * getFitEsts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getFitEsts
### Title: Get Estimates from icenReg Regression Model
### Aliases: getFitEsts

### ** Examples

	simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
                         inspections = 6, inspectLength = 1)
	fit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2,
                  data = simdata)
	new_data <- data.frame(x1 = c(1,2), x2 = c(-1,1))
	rownames(new_data) <- c('grp1', 'grp2')
	
	estQ <- getFitEsts(fit, new_data, p = c(.25, .5, .75))
	
	estP <- getFitEsts(fit, q = 400)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getFitEsts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("getSCurves")
### * getSCurves

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getSCurves
### Title: Get Estimated Survival Curves from Semi-parametric Model for
###   Interval Censored Data
### Aliases: getSCurves

### ** Examples

	set.seed(1)

	sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
	                      shape = 2, scale = 2,
	                      inspections = 6, inspectLength = 1)
	fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data, bs_samples = 0)	

	new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
	#want to fit survival curves with above covariates
	rownames(new_data) <- c('group 1', 'group 2')
	#getSCurves will name the survival curves according to rownames

	curveInfo <- getSCurves(fit, new_data)
	xs <- curveInfo$Tbull_ints
	#Extracting Turnbull intervals
	sCurves <- curveInfo$S_curves
	#Extracting estimated survival curves
	
	plot(xs[,1], sCurves[[1]], xlab = 'time', ylab = 'S(t)', 
	     type = 's', ylim = c(0,1),
	     xlim = range(as.numeric(xs), finite = TRUE))
	#plotting upper survival curve estimate
	lines(xs[,2], sCurves[[1]], type = 's')
	#plotting lower survival curve estimate
	
	lines(xs[,1], sCurves[[2]], col = 'blue', type = 's')
	lines(xs[,2], sCurves[[2]], col = 'blue', type = 's')
	#plotting upper and lower survival curves for group 2
	
	# Actually, all this plotting is a unnecessary: 
	# plot(fit, new_data) will bascially do this all
	# But this is more of a tutorial in case custom
	# plots were desired



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getSCurves", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ic_par")
### * ic_par

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ic_par
### Title: Parametric Regression Models for Interval Censored Data
### Aliases: ic_par

### ** Examples

	data(miceData)
	
	logist_ph_fit <- ic_par(Surv(l, u, type = 'interval2') ~ grp, 
	                        data = miceData, dist = 'loglogistic')

	logist_po_fit <- ic_par(Surv(l, u, type = 'interval2') ~ grp, 
	                        data = miceData, dist = 'loglogistic', model = 'po')

	summary(logist_ph_fit)
	summary(logist_po_fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ic_par", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ic_sp")
### * ic_sp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ic_sp
### Title: Semi-Parametric models for Interval Censored Data
### Aliases: ic_sp vcov.icenReg_fit summary.icenReg_fit plot.icenReg_fit

### ** Examples

	set.seed(1)

	sim_data <- simIC_weib(n = 500, inspections = 5, inspectLength = 1)
	ph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data)	
	# Default fits a Cox-PH model
	
	summary(ph_fit)		
	# Regression estimates close to true 0.5 and -0.5 values


	new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
	rownames(new_data) <- c('group 1', 'group 2')
	plot(ph_fit, new_data)
	# plotting the estimated survival curves

	po_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data,
	                model = 'po')
	# fits a proportional odds model
	
	summary(po_fit)
	
	# Not run: how to set up multiple cores
	# library(doParallel)
	# myCluster <- makeCluster(2, type = 'FORK') 
	# registerDoParallel(myCluster)
	# fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2,
	#              data = sim_data, useMCores = TRUE
	#              bs_samples = 500)
	# stopCluster(myCluster)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ic_sp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("miceData")
### * miceData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: miceData
### Title: Lung Tumor Interval Censored Data from Hoel and Walburg 1972
### Aliases: miceData

### ** Examples

	data(miceData)
	
	coxph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ grp, 
	                   bs_samples = 50,	
	                   data = miceData)
	
	#In practice, more bootstrap samples should be used for inference
	#Keeping it quick for CRAN testing purposes 
	
	summary(coxph_fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("miceData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("optCliq")
### * optCliq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: optCliq
### Title: Computes the MLE for a Binary Mixture Model.
### Aliases: optCliq cliqOptInfo

### ** Examples

testData <- simBVCen()
#simulate bivariate interval censored data

cliqMat <- MLEcens::reduc(testData, cm = TRUE)$cm
#computes the cliqMat associated with data

cliqFit <- optCliq(cliqMat)
#optimizes the component weights for clique matrix

cliqFit 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("optCliq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simBVCen")
### * simBVCen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simBVCen
### Title: Simulates Bivariate Interval Censored Data
### Aliases: simBVCen

### ** Examples

  testData <- simBVCen()
  #simulate bivariate interval censored data
  
 bvcenFit <- ICNPMLE(testData)
 
 bvcenFit



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simBVCen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simIC_weib")
### * simIC_weib

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simIC_weib
### Title: Simulates interval censored data from regression model with a
###   Weibull baseline
### Aliases: simIC_weib

### ** Examples

	set.seed(1)
	sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3, model = 'ph', 
	                       shape = 2, scale = 2, inspections = 6, inspectLength = 1)
	#simulates data from a cox-ph with beta weibull distribution.
	
	diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data, model = 'po')
	diag_covar(Surv(l, u, type = 'interval2') ~ x1 + x2, data = sim_data, model = 'ph')
	#'ph' fit looks better than 'po'; the difference between the transformed survival
	#function looks more constant



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simIC_weib", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
