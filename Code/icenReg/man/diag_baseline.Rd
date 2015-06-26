\name{diag_baseline}
\alias{diag_baseline}
\title{Compare parametric baseline distributions with semi-parametric baseline}
\usage{
diag_baseline(object, data, model = "ph", weights = NULL,
              dists = c("exponential", "weibull", "gamma", 
                        "lnorm", "loglogistic"),
              max_n_use = 10000, cols = NULL, 
              lgdLocation = "bottomleft", useMidCovars = T)}
\arguments{
  \item{object}{Either a formula or a model fit with \code{ic_sp} or \code{ic_par}}
  \item{data}{data. Unnecessary if \code{object} is a fit}
  \item{model}{type of model. Choices are \code{'ph'} or \code{'po'}}
  \item{dists}{parametric baseline fits}	
  \item{max_n_use}{maximum number of observations used to fit semi-parametric model}
  \item{cols}{colors of baseline distributions}
  \item{weights}{case weights}
  \item{lgdLocation}{where legend will be placed. See \code{?legend} for more details}
  \item{useMidCovars}{Should the distribution plotted be for covariates = mean values instead of 0}
}
\description{Creates plots to diagnosis fit of different choices of parametric baseline model. Plots the semi paramtric model (which only uses up to \code{max_n_use} subjects from the dataset for speed) against different choices of parametric models. 
  }
\details{
	If \code{useMidCovars = T}, then the survival curves plotted are for fits with the mean covariate value, rather than 0. This is because often the baseline distribution (i.e. with all covariates = 0) will be far away from the majority of the data.
	
	For larger datasets (n > 10,000), the semi-parametric model can be become slow. Because of this, \code{diag_covar} allows the user to select \code{max_n_use} samples from their data to assess the fit of their models. 

}
\examples{
	 data(essIncData_small)
	 useData <- essIncData_small
	 
	 #Note to user: suggest replacing useData with essIncData
	 #instead of essIncData_small. Using small dataset to quickly 
	 #pass CRAN tests
	 
	 fit_po <- ic_sp(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry,
	                 data = useData, bs_samples = 0, model = 'po')
	
	 diag_baseline(fit_po)

	 #Weibull model appears best fit
}
\author{
  Clifford Anderson-Bergman
}