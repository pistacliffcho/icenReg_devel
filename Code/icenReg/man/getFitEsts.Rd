\name{getFitEsts}
\alias{getFitEsts}
\title{Get Estimates from icenReg Regression Model}
\usage{
  getFitEsts(fit, newdata, p, q) 
}
\arguments{
  \item{fit}{model fit with \code{ic_par} or \code{ic_sp}}
  \item{newdata}{data.frame containing covariates}
  \item{p}{percentiles}
  \item{q}{quantiles}
}
\description{ Gets estimates from a \code{ic_par} or \code{ic_sp} object. Provided estimates conditional on regression parameters found in \code{newdata}.
	 }
 \details{ 	
	If \code{newdata} is left blank, baseline estimates will be returned (i.e. all covariates = 0). If \code{p} is provided, will return the estimated F^{-1}(p | x). If \code{q} is provided, will return the estimated F(q | x). If neither \code{p} nor \code{q} are provided, the estimated conditional median is returned.
	 	
 	For \code{ic_par} fits, it is worth noting that F^{-1}(p | x) is approximated using R's \code{optimize} function, so there may be some numerical error. 
 	
 	In the case of \code{ic_sp}, the MLE of the baseline survival is not necessarily unique, as probability mass is assigned to disjoint Turnbull intervals, but the likelihood function is indifferent to how probability mass is assigned within these intervals. In order to have a well defined estimate returned, we assume probability is assigned uniformly in these intervals. In otherwords, we return *a* maximum likelihood estimate, but don't attempt to characterize *all* maximum likelihood estimates with this function. If that is desired, all the information needed can be extracted with \code{getSCurves}. 
 	 }
\examples{
	simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
                         inspections = 6, inspectLength = 1)
	fit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2,
                  data = simdata)
	new_data <- data.frame(x1 = c(1,2), x2 = c(-1,1))
	rownames(new_data) <- c('grp1', 'grp2')
	
	estQ <- getFitEsts(fit, new_data, p = c(.25, .5, .75))
	
	estP <- getFitEsts(fit, q = 400)
}
\author{
  Clifford Anderson-Bergman
}