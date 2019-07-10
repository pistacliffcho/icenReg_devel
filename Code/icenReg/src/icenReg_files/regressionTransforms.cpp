//
//  regressionTransforms.cpp
//  
//
//  Created by Cliff Anderson Bergman on 11/17/15.
//
//

#include "regressionTransforms.h"

SEXP s_regTrans(SEXP times, SEXP etas,
                SEXP bli, SEXP regType, SEXP baseType){
    condProbCal rtObj(regType, baseType, bli);
    if(rtObj.isOK == false){
        Rprintf("s_regTrans not okay for some reason\n");
        return(R_NilValue);
    }
    int k = LENGTH(times);
    if(k != LENGTH(etas)){
        Rprintf("warning: LENGTH(times) != LEGNTH(etas). Quiting\n");
        return(R_NilValue);
    }
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, k));
    double* c_time = REAL(times);
    double* c_etas = REAL(etas);
    double val;
    for(int i = 0; i < k; i++){
        val = rtObj.getBaseSurv(c_time[i], bli);
        REAL(ans)[i] = 1.0 - rtObj.transformSurv(val, c_etas[i]);
    }
    UNPROTECT(1);
    return(ans);
}


SEXP q_regTrans(SEXP p, SEXP etas,
                SEXP bli, SEXP regType, SEXP baseType){
    
    condProbCal rtObj(regType, baseType, bli);
    if(rtObj.isOK == false){
        Rprintf("s_regTrans not okay for some reason\n");
        return(R_NilValue);
    }
    int k = LENGTH(p);
    if(k != LENGTH(etas)){
        Rprintf("warning: LENGTH(times) != LEGNTH(etas). Quiting\n");
        return(R_NilValue);
    }
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, k));
    double val;
    double* c_p = REAL(p);
    double* c_etas = REAL(etas);
    for(int i = 0; i < k; i++){
        val = rtObj.transform_p(c_p[i], c_etas[i]);
        REAL(ans)[i] = rtObj.getBaseQ(val, bli);
    }
    UNPROTECT(1);
    return(ans);
}


condProbCal::condProbCal(SEXP regType, SEXP baseType, SEXP bli){
    int cRegType = INTEGER(regType)[0];
    isOK = false;
    if(cRegType == 1){
        transformSurv = &propHazTrans;
        transform_p = &transform_p_ph;
    }
    else if(cRegType == 2){
        transformSurv = &propOddsTrans;
        transform_p  = &transform_p_po;
    }
    else if(cRegType == 0){
    	transformSurv = &noTrans;
    	transform_p = &transform_p_none;
    }
    else{
        Rprintf("warning: invalid regType selected. Setting to Cox PH\n");
        transformSurv = &propHazTrans;
        transform_p = &transform_p_ph;
    }
    
    int cBaseType = INTEGER(baseType)[0];
    if(cBaseType == 1){
        getBaseSurv = &getGammaSurv;
        getBaseQ    = &getGammaQ;
        if(LENGTH(bli) == 2){
            isOK = true;
        }
    }
    else if(cBaseType == 2){
        getBaseSurv = &getWeibSurv;
        getBaseQ    = &getWeibQ;
        if(LENGTH(bli) == 2){
            isOK = true;
        }
    }
    else if(cBaseType == 3){
        getBaseSurv = &getLogNormSurv;
        getBaseQ    = &getLogNormQ;
        if(LENGTH(bli) == 2){
            isOK = true;
        }
    }
    else if(cBaseType == 4){
        getBaseSurv = &getExpSurv;
        getBaseQ    = &getExpQ;
        if(LENGTH(bli) == 1){
            isOK = true;
        }
    }
    else if(cBaseType == 5){
        getBaseSurv = &getLgLgsticSurv;
        getBaseQ    = &getLgLgsticQ;
        if(LENGTH(bli) == 2){
            isOK = true;
        }
    }
    else if(cBaseType == 6){
        getBaseSurv = &getGenGammaSurv;
        getBaseQ    = &getGenGammaQ;
        if(LENGTH(bli) == 3){
            isOK = true;
        }
    }
    else if(cBaseType == 0){
        getBaseSurv = &getNonParSurv;
        getBaseQ    = &getNonParQ;
        if(LENGTH(bli) == 2){
            isOK = true;
        }
    }
}





double getNonParSurv(double t, SEXP SC){	
    SEXP tb_ints = VECTOR_ELT(SC, 0);
    SEXP svals = VECTOR_ELT(SC, 1);
    PROTECT(tb_ints);
    PROTECT(svals);
    UNPROTECT(2);
    
    int k = LENGTH(svals);
    if(k != LENGTH(tb_ints)/2){
        Rprintf("LENGTH(tb_ints) = %d, LENGTH(svals) = %d\n", LENGTH(tb_ints), k);
        return(0.0);
    }
    
    double* tb_ptr = REAL(tb_ints);
    double* svals_ptr = REAL(svals);
    
    int ind = 0;
    while(ind < k && t > tb_ptr[ind + k]){ ind++; }
    if(ind == k){ return(0.0); }
    if(ind == 0){ return(1.0); }
    if(t < tb_ptr[ind]){ return(svals_ptr[ind - 1]); }
    double intLength = tb_ptr[ind + k] - tb_ptr[ind];
    double ans = svals_ptr[ind-1];
	
	if(intLength == 0){
		ans = svals_ptr[ind];
		if(ans < 0){
			Rprintf("warning: ans < 0. intLength == 0\n");
			return(ans);
		}
		return(svals_ptr[ind]);
	}
	if(intLength < 0){
		Rprintf("Warning: intLength < 0 when getting survival probabilities\n");
		return(-1.0);
	}		   
    double t_diff = t - tb_ptr[ind];
	if(t_diff < 0){
		Rprintf("warning: t_diff < 0. t = %f, tb_ptr[ind] = %f\n", t, tb_ptr[ind]);
	}
	
	if(t_diff > (tb_ptr[ind + k] - tb_ptr[ind])){ t_diff = intLength;} 
	//t may be outside this turnbull interval but not into the next turnbull interval
    double pStep = svals_ptr[ind-1] - svals_ptr[ind];
	ans -= pStep * t_diff/intLength; 
	
	if(ans < 0 || ans > 1){
		Rprintf("ans < 0 || ans > 1. t_diff = %f, pStep = %f, intLen = %f, ind = %d, k = %d \n",
				t_diff, pStep, intLength, ind, k);
		return(ans);
	}
	
    return(ans);
}

double getNonParQ(double p, SEXP SC){
    double s = 1-p;
    SEXP tb_ints = VECTOR_ELT(SC, 0);
    SEXP svals = VECTOR_ELT(SC, 1);
    PROTECT(tb_ints);
    PROTECT(svals);
    UNPROTECT(2);
    
    int k = LENGTH(svals);
    if(k != LENGTH(tb_ints)/2){
        Rprintf("LENGTH(tb_ints) = %d, LENGTH(svals) = %d\n", LENGTH(tb_ints), k);
        return(0.0);
    }
    double* tb_ptr = REAL(tb_ints);
    double* svals_ptr = REAL(svals);

    int ind = 0;
    while( svals_ptr[ind] > s && ind < k){ ind++; }
    if(ind == 0){ return(tb_ptr[0]); }
    if(ind == k){ return(tb_ptr[ind + k - 1]); }

    double intLength = tb_ptr[ind + k] - tb_ptr[ind];
    double s_diff = svals_ptr[ind-1] - s;
    double pLength = svals_ptr[ind-1] - svals_ptr[ind];
    
    double ans = tb_ptr[ind] + intLength * s_diff/pLength;

    return(ans);
}
