\name{diag_covar}
\alias{diag_covar}
\title{Evaluate covariate effect for regression model}
\usage{
diag_covar(object, varName, 
           data, model, weights = NULL,
           yType = 'meanRemovedTransform',
           factorSplit = TRUE, 
           numericCuts, col, 
           xlab, ylab, main, 
           lgdLocation = NULL) }
\arguments{
  \item{object}{Either a formula or a model fit with \code{ic_sp} or \code{ic_par}}
  \item{varName}{covariate to split data on. If left blank, will split on each covariate}
  \item{data}{data. Unnecessary if \code{object} is a fit}
  \item{model}{type of model. Choices are \code{'ph'} or \code{'po'}}
  \item{weights}{case weights}
  \item{yType}{type of plot created. See details}
  \item{factorSplit}{Should covariate be split as a factor (i.e. by levels)}
  \item{numericCuts}{If \code{fractorSplit == FALSE}, cut points of covariate to stratify data on} 
  \item{col}{colors of each subgroup plot. If left blank, will auto pick colors}
  \item{xlab}{label of x axis}
  \item{ylab}{label of y axis}
  \item{main}{title of plot}
  \item{lgdLocation}{where legend should be placed. See details}
}
\description{Creates plots to diagnosis fit of covariate effect in a regression model. For a given variable, stratifies the data across different levels of the variable and adjusts for all the other covariates included in \code{fit} and then plots a given function to help diagnosis where covariate effect follows model assumption (i.e. either proportional hazards or proportional odds). See \code{details} for descriptions of the plots. 
	
If \code{varName} is not provided, will attempt to figure out how to divide up each covariate and plot all of them, although this may fail. 
  }
\details{
	For the Cox-PH and proportional odds models, there exists a transformation of survival curves such that the difference should be constant for subjects with different covariates. In the case of the Cox-PH, this is the log(-log(S(t|X))) transformation, for the proporitonal odds, this is the log(S(t|X) / (1 - S(t|X))) transformation. 
	
	The function diag_covar allows the user to easily use these transformations to diagnosis whether such a model is appropriate. In particular, it takes a single covariate and stratifies the data on different levels of that covariate. Then, it fits the semi-parametric regression model (adjusting for all other covariates in the data set) on each of these stratas and extracts the baseline survival function. If the stratified covariate does follow the regression assumption, the difference between these transformed baseline survival functions should be approximately constant. 
	
	To help diagnosis, the default function plotted is the transformed survival functions, with the overall means subtracted off. If the assumption holds true, then the difference between the plotted lines should be approximately constant. Other choices of \code{yType}, the function to plot, are \code{"transform"}, which is the transformed functions without the means subtracted and \code{"survival"}, which is the baseline survival distribution is plotted for each strata. 
	
	Currently does not support stratifying covariates that are involved in an interaction term. 
	
	For variables that are factors, it will create a strata for each level of the covariate, up to 20 levels. If \code{factorSplit == FALSE}, will divide up numeric covariates according to the cuts provided to numericCuts. 
	
	\code{lgdLocation} is an argument placed to \code{legend} dictating where the legend will be placed. If \code{lgdLocation = NULL}, will use standard placement given \code{yType}. See \code{?legend} for more details. 
}
\examples{
	 data(essIncData_small)
	 useData <- essIncData_small
	 
	 #Note to user: suggest replacing useData with essIncData
	 #instead of essIncData_small. Using small dataset to quickly 
	 #pass CRAN tests
	 par(mfrow = c(1,2))
	 	 
	 diag_covar(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 
		         data = useData, model = 'po')
		
	 diag_covar(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 
		         data = useData, model = 'ph')}
\author{
  Clifford Anderson-Bergman
}