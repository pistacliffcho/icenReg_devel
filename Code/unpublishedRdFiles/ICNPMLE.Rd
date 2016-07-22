\name{ICNPMLE}
\alias{ICNPMLE}
\alias{IC_NPMLE}
\docType{package}

\title{Computes the NPMLE for Univariate or Bivariate Interval Censored Data}

\author{Clifford Anderson-Bergman

Also uses Marloes Maathius's \code{MLEcens::reduc} function for calculation of the clique matrix.}

\description{
Computes the MLE for a Interval Censored Data with a squeezing EM algorithm (\emph{much} faster than the standard EM). Accepts either univariate interval censored data (where \code{times} is an n x 2 matrix with \code{times[,1]} being the left side of the interval and \code{times[,2]} is the right side), or bivariate interval censored data (where \code{times} is an n x 4 matrix with \code{times[,1:2]} being left and right side of the interval for event 1 and \code{times[,3:4]} being the left and right side of the interval for event 2).
}
\usage{
 ICNPMLE(times, B = c(1,1), max.inner = 100, max.outer = 100, tol = 1e-10)
}
\arguments{
  \item{times}{either an  n x 2 or n x 4 \code{data.frame} or \code{matrix} of censoring intervals}
  \item{B}{A vector indicating whether each end of the intervals are open (0) or closed (1). Alternatively, this could be an n x 2 or n x 4 matrix of indicators for each individual interval}
  \item{max.inner}{number of inner loops used in optimization algorithm}
  \item{max.outer}{number of outer loops used in optimization algorithm}
  \item{tol}{numerical tolerance}
}

\references{
Anderson-Bergman, C., (2014) Semi- and non-parametric methods for interval censored data with shape constraints, Ph.D. Thesis

Yu, Y., (2010), Improved EM for Mixture Proportions with Applications to Nonparametric ML Estimation for Censored Data, \emph{preprint}

Maathuis, M., (2005). Reduction algorithm for the NPMLE for the distribution function of bivariate interval censored data. \emph{Journal of Computational and Graphical Statistics} 
Vol 14 pp 252\- 262}
\examples{
simData <- simBVCen(500)

fit <- ICNPMLE(simData)

fit
}
