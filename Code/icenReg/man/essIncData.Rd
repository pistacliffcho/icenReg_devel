\name{essIncData}
\alias{essIncData}
\title{Interval Censored Income Data from European Social Survey}
\usage{data(essIncData)}
\format{A data frame with 6712 rows and 4 variables	
	\itemize{
	  \item{\code{cntry}}{   Country}
	  \item{\code{eduLevel}}{   Categorical variable for number of years of reported education}
	  \item{\code{inc_l}}{   Lower limit of reported income level (in Euros)}
	  \item{\code{inc_u}}{   Upper limit of reported income level (in Euros)}
	  }
}
\description{ 
	Dataset containing a sample from the European Social Survey.
	
	In the European Social Survey, income is reported up a to discrete intervals. Within each country, these intervals are disjoint (i.e. [0, 1k), [1k,2k)...). However, across countries, the intervals are not disjoint and so interval censored methods should be used to compare subjects across different countries.
 	}
\examples{
	data(essIncData)
	
	lnormFit <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, 	
	                   data = essIncData,
	                   model = 'po',
	                   dist = 'loglogistic')
	
	summary(lnormFit)
}
\source{
	ESS Round 5: European Social Survey Round 5 Data (2010). Data file edition 3.2. Norwegian Social Science Data Services, Norway\- Data Archive and distributor of ESS data.
}