//
//  icenReg.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

//#include "Eigen_local/Dense"
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>
#include <Rcpp.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		

using namespace std;

#include "icenReg_files/basicUtilities.cpp"
#include "icenReg_files/ic_par.h"
#include "icenReg_files/ic_par_aft.h"
#include "icenReg_files/ic_par.cpp"
#include "icenReg_files/ic_sp_ch.cpp"
#include "icenReg_files/bivariateNPMLE.cpp"
#include "icenReg_files/ic_sp_gradDescent.cpp"

#include "icenReg_files/experimentalCode.cpp"
#include "icenReg_files/regressionTransforms.cpp"
#include "icenReg_files/regressionTransforms_2.cpp"

#include "icenReg_files/EMICM.h"
#include "icenReg_files/EMICM.cpp"
#include "icenReg_files/ic_par_aft.cpp"

#include "icenReg_files/ic_bayes.cpp"

/*
// [[Rcpp::export]]
Rcpp::List ic_par(SEXP R_s_t,     // this is a vector of times for which the baseline survival will be calculated
                SEXP R_d_t,     // this is a vector of times for which the baseline density will be calculated
                SEXP covars,    // matrix of covariates. Will need to be ordered according to group
                SEXP uncenInd,  // a n_1 x 2 matrix of indicators. First will be index to look up corresponding d
                                // (i.e. for d_t), second is index to look up s_t
                SEXP gicInd,    // a n_2 x 2 matrix of indicators. First will be index to look up left side of interval,
                                // second is to look up right side of interval
                SEXP lInd,      // n_3 vector of indicators of right side of interval
                SEXP rInd,      // n_4 vector of indicators of left side of interval
                SEXP parType,   // integer of paramteric type. 1 = gamma, 2 = weib, 3 = lnorm, 4 = exp, 5 = loglogistic
                SEXP linkType,  // integer of link type. 1 = proportional odds, 2 = proportional hazards
                SEXP outHessian, // hessian matrix at MLE. Easier to pass this than to create it in C++
                SEXP R_w         //
                );
*/

// [[Rcpp::export]]
Rcpp::List ic_parList(Rcpp::List R_parList);

// [[Rcpp::export]]
Rcpp::List R_ic_bayes(Rcpp::List R_bayesList, Rcpp::Function priorFxn, 
					  Rcpp::List R_ic_parList);
					  
//[[Rcpp::export]]
Rcpp::NumericVector computeConditional_p(Rcpp::NumericVector q,
										 Rcpp::NumericVector etas,
										 Rcpp::NumericMatrix baselineParams,
										 Rcpp::CharacterVector reg_model,
										 Rcpp::CharacterVector base_dist);
										 
//[[Rcpp::export]]
Rcpp::NumericVector computeConditional_q(Rcpp::NumericVector p,
										 Rcpp::NumericVector etas,
										 Rcpp::NumericMatrix baselineParams,
										 Rcpp::CharacterVector reg_model,
										 Rcpp::CharacterVector base_dist);
