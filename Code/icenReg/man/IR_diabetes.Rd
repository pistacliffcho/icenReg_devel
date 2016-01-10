\name{IR_diabetes}
\alias{IR_diabetes}
\title{Interval censored time from diabetes onset to diabetic nephronpathy}
\usage{data(IR_diabetes)}
\format{A data frame with 731 rows and 3 variables	
  \itemize{
    \item{\code{left}}{   left side of observation interval}
    \item{\code{right}}{   right side of observation interval}
    \item{\code{gender}}{   gender of subject}
  }
}
\description{
  Data set contains interval censored survival time for time from onset of
  diabetes to to diabetic nephronpathy. Identical to the \code{diabetes}
  dataset found in the package \code{glrt}. 
  }
\examples{
data(IR_diabetes)
fit <- ic_par(cbind(left, right) ~ gender, 
              data = IR_diabetes,
              model = "po",
              dist = "loglogistic")
}
\source{
  Borch-Johnsens, K, Andersen, P and Decker, T (1985).
  "The effect of proteinuria on relative mortality in Type I (insulin-dependent) diabetes mellitus."
  Diabetologia, 28, 590-596.
}