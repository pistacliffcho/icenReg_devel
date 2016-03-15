# icenReg_devel
Repository for development of R package icenReg

This is an R-package for regression models for interval censored data. 

Currently supports semi-parametric proportional hazards or odds models,
along with fully parametric versions of the same models. Does not fit an AFT model, as a
parametric model of this sort can be fit using survival's survreg() function.
It can also fit the non-parametric maximum likelihood estimate (NPMLE) for interval censored data.

Also contains tools for diagnostic checks of models and imputation of censored data.

If you use the package and find a bug, want a feature added or can think of someway to
improve it, please let me know!
