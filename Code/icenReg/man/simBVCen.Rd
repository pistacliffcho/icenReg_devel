\name{simBVCen}
\alias{simBVCen}
\docType{package}
\title{
  Simulates Bivariate Interval Censored Data
}
\description{
  Simulates Bivariate Interval Censored Data
}
\usage{
  simBVCen(n = 1000)
}
\arguments{
  \item{n}{number of observations simulated}

}

\examples{
  testData <- simBVCen()
  #simulate bivariate interval censored data
  
 bvcenFit <- bivariateNPMLE(testData)
 
 bvcenFit
}
