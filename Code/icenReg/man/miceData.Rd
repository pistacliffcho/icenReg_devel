\name{miceData}
\alias{miceData}
\title{Lung Tumor Interval Censored Data from Hoel and Walburg 1972}
\usage{data(miceData)}
\format{A data frame with 144 rows and 3 variables	
	\itemize{
	  \item{\code{l}}{   left side of observation interval}
	  \item{\code{u}}{   right side of observation interval}
	  \item{\code{grp}}{   Group for mouse. Either ce (conventional environment) or ge (grem-free environment)}
	  }
}
\description{
 	RFM mice were sacrificed and examined for lung tumors. This resulted in current status interval censored data: if the tumor was present, this implied left censoring and if no tumor was present this implied right censoring. Mice were placed in two different groups: conventional environment or germ free environment.}
\examples{
data(miceData)
	
coxph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ grp, 
                   bs_samples = 50,	
                   data = miceData)
	
#In practice, more bootstrap samples should be used for inference
#Keeping it quick for CRAN testing purposes 
	
summary(coxph_fit)
}
\source{
	Hoel D. and Walburg, H.,(1972), Statistical analysis of survival experiments, \emph{The Annals of Statistics}, 18, 1259-1294 
}