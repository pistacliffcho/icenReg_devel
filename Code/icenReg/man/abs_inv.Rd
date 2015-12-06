\name{abs_inv}
\alias{abs_inv}
\title{Loss Function for interval censored cross validation}
\usage{
  abs_inv(pred, t_val)
}
\arguments{
  \item{pred}{predicted value}
  \item{t_val}{true value}
}
\description{
  A loss function for survival data. Equal to 
  
  \code{mean(abs(1/(pred+1) - 1/(t_val+1 ) ) ) }. 
  
  This function heavily penalizes inaccurate predictions that are closer to the 
  origin than those that are farther way. In otherwords, it will favor a model
  that accurately predicts high risk subjects over a model that accurately predicts
  loss risk subjects. 
}

\author{
  Clifford Anderson-Bergman
}