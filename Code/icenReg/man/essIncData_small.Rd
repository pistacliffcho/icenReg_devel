\name{essIncData_small}
\alias{essIncData_small}
\title{Interval Censored Income Data from European Social Survey}
\usage{data(essIncData_small)}
\format{A data frame with 500 rows and 4 variables	
	\itemize{
	  \item{\code{cntry}}{   Country}
	  \item{\code{eduLevel}}{   Categorical variable for number of years of reported education}
	  \item{\code{inc_l}}{   Lower limit of reported income level (in Euros)}
	  \item{\code{inc_u}}{   Upper limit of reported income level (in Euros)}
	  }
}
\description{ 
	Dataset containing a sample from the European Social Survey. This is a small subsample from the dataset \code{essIncData}, used only to run the example very quickly for CRAN. Using the full dataset, running the examples should still be fairly fast; at most a few seconds. 
	
	In the European Social Survey, income is reported up a to discrete intervals. Within each country, these intervals are disjoint (i.e. [0, 1k), [1k,2k)...). However, across countries, the intervals are not disjoint and so interval censored methods should be used to compare subjects across different countries.
 	}
\examples{
	data(essIncData_small)
	
	lnormFit <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 	
	                   data = essIncData_small,
	                   model = 'po',
	                   dist = 'loglogistic')
	
	summary(lnormFit)
}
\source{
	ESS Round 5: European Social Survey Round 5 Data (2010). Data file edition 3.2. Norwegian Social Science Data Services, Norway\- Data Archive and distributor of ESS data.
}