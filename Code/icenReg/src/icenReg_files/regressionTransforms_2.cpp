//
//  regressionTransforms_2.cpp
//  
//
//  Created by Cliff Anderson Bergman on 11/17/15.
//
//

#include "regressionTransforms_2.h"


Rcpp::NumericVector computeConditional_p(Rcpp::NumericVector q,
										 Rcpp::NumericVector etas,
										 Rcpp::NumericMatrix baselineParams,
										 Rcpp::CharacterVector reg_model,
										 Rcpp::CharacterVector base_dist){
	condProbCal_2 calObj(reg_model, base_dist);
	int nRow = baselineParams.nrow();
	Rcpp::NumericVector ans(nRow);
	std::vector<double> these_baselines;
	for(int i = 0; i < nRow; i++){
		these_baselines = getRow(i, baselineParams);
		ans[i] = calObj.calc_p(q[i], etas[i], these_baselines);
	}
	return(ans);
}

Rcpp::NumericVector computeConditional_q(Rcpp::NumericVector p,
										 Rcpp::NumericVector etas,
										 Rcpp::NumericMatrix baselineParams,
										 Rcpp::CharacterVector reg_model,
										 Rcpp::CharacterVector base_dist){
	condProbCal_2 calObj(reg_model, base_dist);
	int nRow = baselineParams.nrow();
	Rcpp::NumericVector ans(nRow);
	std::vector<double> these_baselines;
	for(int i = 0; i < nRow; i++){
		these_baselines = getRow(i, baselineParams);
		ans[i] = calObj.calc_q(p[i], etas[i], these_baselines);
	}
	return(ans);
}

double condProbCal_2::calc_p(double q, double nu, std::vector<double> &bli){
	if(isAFT){ q  = q/nu; }
	double base_s = getBase_s(q, bli);
	double ans    = baseSurv_2_condSurv(base_s, nu);
	ans = 1.0 - ans;
	return(ans);
}

double condProbCal_2::calc_q(double p, double nu, std::vector<double> &bli){
	double s = 1.0 - p;
	double base_s = condSurv_2_baseSurv(s, nu);
	double ans    = getBase_q(base_s, bli);
	if(isAFT){ ans *= nu; }
	return(ans);
}


condProbCal_2::condProbCal_2(Rcpp::CharacterVector regType, 
							 Rcpp::CharacterVector baseType){
    isOK = false;
    isAFT = false;
    if(regType[0] == "ph"){
    	baseSurv_2_condSurv = &baseSurv_2_condSurv_ph;
    	condSurv_2_baseSurv = &condSurv_2_baseSurv_ph;
    }
    else if(regType[0] == "po"){
    	baseSurv_2_condSurv = &baseSurv_2_condSurv_po;
    	condSurv_2_baseSurv = &condSurv_2_baseSurv_po;
    }
    else if(regType[0] == "aft"){
    	baseSurv_2_condSurv = &noTrans;
    	condSurv_2_baseSurv = &noTrans;
    	isAFT = true;
    }
    else{
    	Rcpp::stop("regType not recongized");
/*        Rprintf("warning: invalid regType selected. Setting to Cox PH\n");
    	baseSurv_2_condSurv = &baseSurv_2_condSurv_ph;
    	condSurv_2_baseSurv = &condSurv_2_baseSurv_ph; */
    }
    
    if(baseType[0] == "gamma"){
        getBase_s   = &getGammaSurv;
        getBase_q   = &getGammaQ;
    }
    else if(baseType[0] == "weibull"){
        getBase_s   = &getWeibSurv;
        getBase_q   = &getWeibQ;
    }
    else if(baseType[0] == "lnorm"){
        getBase_s    = &getLogNormSurv;
        getBase_q    = &getLogNormQ;
    }
    else if(baseType[0] == "exponential"){
        getBase_s    = &getExpSurv;
        getBase_q    = &getExpQ;
    }
    else if(baseType[0] == "loglogistic"){
        getBase_s    = &getLgLgsticSurv;
        getBase_q    = &getLgLgsticQ;
    }
    else if(baseType[0] == "generalgamma"){
        getBase_s    = &getGenGammaSurv;
        getBase_q    = &getGenGammaQ;
    }
    else if(baseType[0] == "np"){
        Rcpp::stop("conProbCal_2 currently does not support non/semi-parametric models");
    }
    else{
    	Rcpp::stop("baseType not recongized");
    }
}
