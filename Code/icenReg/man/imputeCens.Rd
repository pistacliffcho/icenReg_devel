\name{imputeCens}
\alias{imputeCens}
\title{Impute Interval Censored Data from icenReg Regression Model}
\usage{
imputeCens(fit, newdata = NULL, 
           imputeType = "fullSample",
          numImputes = 5)
}
\arguments{
  \item{fit}{model fit with \code{ic_par} or \code{ic_sp}}
  \item{newdata}{data.frame containing covariates and censored intervals. If blank, will use data from model}
  \item{imputeType}{type of imputation. See details for options}
  \item{numImputes}{Number of imputations (ignored if \code{imputeType = "median"}) }
}
\description{ 
  Imputes censored responses from data. 
}
\details{ 	
If \code{newdata} is left blank, will provide estimates for original data set. 
  
There are several options for how to impute. \code{imputeType = 'median'} 
imputes the median time, conditional on the response interval, covariates and 
regression parameters at the MLE. To get random imputations without accounting
for error in the estimated parameters \code{imputeType ='fixedParSample'} takes a 
random sample of the response variable, conditional on the response interval, 
covariates and estimated parameters at the MLE. Finally, 
\code{imputeType = 'fullSample'} first takes a random sample of the coefficients,
(assuming asymptotic normality) and then takes a random sample 
of the response variable, conditional on the response interval, 
covariates, and the random sample of the coefficients. 


 }
\examples{

simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
                      inspections = 6, inspectLength = 1)
  
fit <- ic_par(cbind(l, u) ~ x1 + x2,
                   data = simdata)

imputedValues <- imputeCens(fit)
}
\author{
  Clifford Anderson-Bergman
}