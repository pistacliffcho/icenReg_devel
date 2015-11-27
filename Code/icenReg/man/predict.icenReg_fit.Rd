\name{predict.icenReg_fit}
\alias{predict.icenReg_fit}
\alias{predict}
\title{Predictions from icenReg Regression Model}
\usage{
\method{predict}{icenReg_fit}(object, type = "response", 
          newdata = NULL, ...) 
}
\arguments{
  \item{object}{model fit with \code{ic_par} or \code{ic_sp}}
  \item{type}{type of prediction. Options include "lp", "response"}
  \item{newdata}{data.frame containing covariates}
  \item{...}{other arguments (will be ignored)}
}
\description{ 
  Gets various estimates from a \code{ic_par} or \code{ic_sp} object. 
	 }
 \details{ 	
	If \code{newdata} is left blank, will provide estimates for original data set. 
	
	For the argument \code{type}, there are two options. \code{"lp"} provides the linear predictor for each subject (i.e. in a Cox-PH, this is the log-hazards ratio, in proportional odds, the log proporitonal odds), \code{"response"} provides the median response value for each subject, *conditional on that subject's covariates, but ignoring their actual response interval*. Use \code{imputeCens} to impute the censored values. 
}
\examples{

simdata <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
                      inspections = 6, inspectLength = 1)
  
fit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2,
                   data = simdata)

imputedValues <- predict(fit)
}
\author{
  Clifford Anderson-Bergman
}