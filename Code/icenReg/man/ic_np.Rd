\name{ic_np}
\alias{ic_np}
\title{Non-Parametric Estimator for Interval Censored Data}
\usage{
  ic_np(formula = NULL, data, maxIter = 1000, tol = 10^-10) 
}
\arguments{
  \item{formula}{formula for stratification. If only one group, can be left blank and data must be entered as n x 2 matrix.}
  \item{data}{a n x 2 matrix or dataframe. See details.}
  \item{maxIter}{maximum iterations}
  \item{tol}{numeric tolerance}
}
\description{
  Fits the non-parametric maximum likelihood estimator (NPMLE) for univariate interval censored data. This is a generalization of the Kaplan-Meier curves that allows for 
interval censoring. Also referred to as the Turnbull estimator.
}

\details{
  \code{data} must be an n x 2 matrix or data.frame containing two columns of data representing left and right sides of the censoring interval, denoted L and R. This allows for left censored (L == 0), right censored (R == inf), uncensored (L == R) along with general interval censored observations. 
  
  The NPMLE is fit using an efficient implementation of the EMICM algorithm. 
} 

\references{
  Turnbull, B. (1976) The empricial distribution with arbitrarily grouped and censored data, \emph{Journal of the Royal Statistical Society B}, vol 38 p290-295

  Wellner, J. A., and Zhan, Y. (1997) A hybrid algorithm for computation of the maximum likelihood estimator from censored data, \emph{Journal of the  American Statistical Association}, Vol 92, pp945-959
  
  Anderson-Bergman, C. (preprint) An efficient implementation of the EMICM algorithm for the interval censored NPMLE
}

\examples{
  data(IR_diabetes)

  fit <- ic_np(cbind(left, right) ~ gender,
               data = IR_diabetes)  

  # Stratifies by gender
  # If only one group, okay to enter data as n x 2 matrix
  
  plot(fit, col = c('blue', 'orange') )
}
\author{
  Clifford Anderson-Bergman
}