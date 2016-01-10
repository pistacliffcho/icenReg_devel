\name{icenReg_cv}
\alias{icenReg_cv}
\title{Get Cross-Validation Statistics from icenReg Regression model}
\usage{
  icenReg_cv(fit, loss_fun = abs_inv, folds = 10, 
             numImputes = 100, useMCore = F)
}
\arguments{
  \item{fit}{model fit with \code{ic_par} or \code{ic_sp}. See details }
  \item{loss_fun}{Loss function. See details}
  \item{folds}{Number of cross validation folds}
  \item{numImputes}{number of imputations of missing responses}
  \item{useMCore}{Should multiple cores be used? See details}
}
\description{
  Computes the cross validation error icenReg regression models by using multiple imputations to 
  get estimates of the censored values and taking the average loss across all imputations. Allows 
  user to provide their own loss function. 
}
\details{
  For interval censored data, the loss function cannot be directly calculated. This is because 
  the response variable is only observed up to an interval and so standard loss functions 
  are not well defined (i.e. (10 - [7, 13])^2 is not a single number. Here the predicted value
  is 10, and the response is known to be between 7 and 13). To handle this, \code{icenReg_cv}
  uses multiple imputations; it first takes a posterior sample of the regression coefficients
  and conditional on these coefficients and the censored interval for each value, it draws a sample
  for each of the censored values. This process is repeated \code{numImputes} times to account
  for the uncertainty in the imputation. 
  
  The cross validation error will be computed using the user provided \code{loss_fun}. This function
  must take two vectors of numeric values; the first being the predicted values and the second being
  the true values. The default function is \code{abs_inv}, which we believe is well suited in most 
  survival analysis settings. It should be noted that the predicted values used in \code{icenReg_cv}
  are the estimated median for a given set of covariates. In many cases, this may not necessarily be
  the prediction that minimizes the loss function, but it is often still a reasonable estimate.
  
  The object \code{fit} must be an icenReg regression model, i.e. either a parametric model 
  (via \code{ic_par}) or a semiparametric model (via \code{ic_sp}). We caution that the 
  semi-parametric models appear to have heavier bias in estimating the out of sample bias, at 
  least with the \code{abs_inv} loss function. The reason for this is that the semi-parametric
  models \code{only} assign probability mass to to \emph{Turnbull intervals}, which in most cases
  will be strictly greater than 0, while most parametric models have possible probability on the real
  line. This biases both the imputations and the estimates of the semi-parametric model to not 
  be close to 0, which is heavily pelanized by the \code{abs_inv} function. By both imputing no 
  values close to 0 and predicting no values close to 0, the \code{ic_sp} model is being overly
  optimistic. In short, caution
  should be used when using the CV criteria to compare parametric and semi-parametric models. 
  
  Multiple cores can be used to speed up the cross validation process. This must be done with 
  the \code{doParallel} package. 
}
\examples{
  ##Not run: somewhat computationally intense
  # should not take more than a few minutes in worst case scenario
  
  #library(splines)
  #library(doParallel)
  #simData <- simIC_weib(n = 1000)
  
  # fit1 <- ic_par(cbind(l,u) ~ x1 + x2, data = simData)  
       # True model
  
  # fit2 <- ic_par(cbind(l,u) ~ x1 + x2, model = 'po', data = simData) 
       # Model with wrong link function
  
  # fit3 <- ic_par(cbind(l,u) ~ bs(x1, df = 6) + x2, data = simData)  
       # Overfit model
  
  # myCluster <- makeCluster(5)
  # registerDoParellel(myCluster)
        #Setting up to use multiple cores
  
  # icenReg_cv(fit1, useMCore = T)
  # icenReg_cv(fit2, useMCore = T)
  # icenReg_cv(fit3, useMCore = T)
  
  # stopCluster(myCluster)
}
\author{
  Clifford Anderson-Bergman
}