\name{ic_par}
\alias{ic_par}
\title{Parametric Regression  Models for Interval Censored Data}
\usage{
  ic_par(formula, data, model = 'ph', dist = 'weibull', weights = NULL) 
}
\arguments{
  \item{formula}{regression formula. Response must be a \code{Surv} object of type \code{'interval2'}}
  \item{data}{dataset}
  \item{model}{What type of model to fit. Current choices are "\code{ph}" (Cox PH) or "\code{po}" (proportional odds)}
  \item{dist}{What baseline parametric distribution to use. See details for current choices}
  \item{weights}{vector of case weights. Not standardized; see details}
}
\description{
 	Fits a parametric regression model for interval censored data. Can fit either a Cox-PH model or a proportional odds model.  
 	
 }
 
\details{
	Currently supported distributions choices are "exponential", "weibull", "gamma", "lnorm" and "loglogistic". 
	
	Does not allow uncensored data points at t = 0 (i.e. where left side of interval == right side == 0), as this will lead to a degenerate estimator for most parametric families. Unlike the current implementation of \code{survreg}, does allow left side of intervals of positive length to 0 and right side to be \code{Inf}. 
	
	In regards to weights, they are not standardized. This means that if weight[i] = 2, this is the equivalent to having two observations with the same values as subject i. 
	
	For numeric stability, if abs(right - left) < 10^-6, observation is considered uncensored rather than interval censored with an extremely small interval. 
} 
  
\examples{
	data(mice_data)
	
	logist_ph_fit <- ic_par(Surv(l, u, type = 'interval2') ~ grp, 
	                        data = mice_data, dist = 'loglogistic')

	logist_po_fit <- ic_par(Surv(l, u, type = 'interval2') ~ grp, 
	                        data = mice_data, dist = 'loglogistic', model = 'po')

	summary(logist_ph_fit)
	summary(logist_po_fit)
}
\author{
  Clifford Anderson-Bergman
}