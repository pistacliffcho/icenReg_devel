\name{ic_sp}
\alias{ic_sp}
\alias{vcov.icenReg_fit}
\alias{summary.icenReg_fit}
\alias{plot.icenReg_fit}
\title{Semi-Parametric models for Interval Censored Data}
\usage{
  ic_sp(formula, data, model = 'ph', weights = NULL,
        bs_samples = 0, useMCores = F, seed = NULL) 
}
\arguments{
  \item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}}
  \item{data}{dataset}
  \item{model}{What type of model to fit. Current choices are "\code{ph}" (Cox PH) or "\code{po}" (proportional odds)}
   \item{weights}{vector of case weights. Not standardized; see details}
  \item{bs_samples}{Number of bootstrap samples used for estimation of standard errors}
  \item{useMCores}{Should multiple cores be used for bootstrap sample? Does not register cluster (see example)}
  \item{seed}{seed for bootstrap. If \code{seed == NULL}, a random seed is still used. See details}
}
\description{
 	Fits a semi-parametric model for interval censored data. Can fit either a Cox-PH model or a proportional odds model.  
 	
 	The covariance matrix for the regression coefficients is estimated via bootstrapping. For large datasets, this can become slow so parallel processing can be used to take advantage of multiple cores via the \code{foreach} package. 
  
 }
 
\details{
	In regards to weights, they are not standardized. This means that if weight[i] = 2, this is the equivalent to having two observations with the same values as subject i. 

	It is very important to note that a random seed is *always* set if bs_samples > 0 (via \code{set.seed(seed)}), which can create problems in simulation studies if the same seed is set in every call to ic_sp during a simulation study. If \code{seed == NULL}, then the starting seed will be \code{round(runif(0, max = 10^8))}, which should be approximately equivalent to not setting a seed. 

  Likelihood function is maximized using an active set algorithm. 
} 
  
\examples{
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
	#              data = sim_data, useMCores = TRUE)	
}
\author{
  Clifford Anderson-Bergman
}