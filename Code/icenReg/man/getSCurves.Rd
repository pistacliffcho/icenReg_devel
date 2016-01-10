\name{getSCurves}
\alias{getSCurves}
\title{Get Estimated Survival Curves from Semi-parametric Model for Interval Censored Data}
\usage{
  getSCurves(fit, newdata) 
}
\arguments{
  \item{fit}{model fit with \code{ic_sp} }
  \item{newdata}{data.frame containing covariates for which the survival curve will be fit to. Rownames from \code{newdata} will be used to name survival curve. If left blank, baseline covariates will be used}
}
\description{
 	Extracts the estimated survival curve(s) from a ic_sp model for interval censored data. Output will be a list with two elements: the first item will be \code{$Tbull_ints}, which is the Turnbull intervals. This is a k x 2 matrix, with the first column being the beginning of the Turnbull interval and the second being the end. This is necessary due to the \emph{representational non-uniqueness}; any survival curve that lies between the survival curves created from the upper and lower limits of the Turnbull intervals will have equal likelihood. See example for proper display of this. The second item is \code{$S_curves}, or the estimated survival probability at each Turnbull interval for individuals with the covariates provided in \code{newdata}. Note that multiple rows may be provided to newdata, which will result in multiple S_curves. 
  }
\examples{
set.seed(1)

sim_data <- simIC_weib(n = 500, b1 = .3, b2 = -.3,
                       shape = 2, scale = 2,
                       inspections = 6, inspectLength = 1)
fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
             data = sim_data, bs_samples = 0)	

new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
#want to fit survival curves with above covariates
	
rownames(new_data) <- c('group 1', 'group 2')
#getSCurves will name the survival curves according to rownames

curveInfo <- getSCurves(fit, new_data)
xs <- curveInfo$Tbull_ints
#Extracting Turnbull intervals
sCurves <- curveInfo$S_curves
#Extracting estimated survival curves
	
plot(xs[,1], sCurves[[1]], xlab = 'time', ylab = 'S(t)', 
     type = 's', ylim = c(0,1),
     xlim = range(as.numeric(xs), finite = TRUE))
#plotting upper survival curve estimate

lines(xs[,2], sCurves[[1]], type = 's')
#plotting lower survival curve estimate
	
lines(xs[,1], sCurves[[2]], col = 'blue', type = 's')
lines(xs[,2], sCurves[[2]], col = 'blue', type = 's')
#plotting upper and lower survival curves for group 2
	
# Actually, all this plotting is a unnecessary: 
# plot(fit, new_data) will bascially do this all
# But this is more of a tutorial in case custom
# plots were desired
}
\author{
  Clifford Anderson-Bergman
}