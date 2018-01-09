//
//  ic_par.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#include "ic_par.h"

double IC_parOpt::calcLike_baseReady(){
    double ans = 0;
    int w_ind = -1;
    int thisSize = uc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_d(d_v[uc[i].d], s_v[uc[i].s], expEta[uc[i].nu])) * w[w_ind] ;
    }
    thisSize = gic.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_s(s_v[gic[i].l], expEta[gic[i].nu])
                   -lnkFn->con_s(s_v[gic[i].r], expEta[gic[i].nu]) ) * w[w_ind];
    }
    thisSize = lc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(1.0 - lnkFn->con_s(s_v[lc[i].r], expEta[lc[i].nu])) * w[w_ind];
    }
    thisSize = rc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(lnkFn->con_s(s_v[rc[i].l], expEta[rc[i].nu])) * w[w_ind];
    }
    
    if(ISNAN(ans)) ans = R_NegInf;
    return(ans);
}


void IC_parOpt::update_dobs_detas(){
    double con0, con_l, con_h, thisEta, thisExpEta, this_d, this_sl, this_sr;
    int w_ind = -1;
    int thisSize = uc.size();
    
    double this_h = h ;
    
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta  = eta[uc[i].nu];
        this_d = d_v[uc[i].d];
        this_sl = s_v[uc[i].s];
        con0 = log(lnkFn->con_d(this_d, this_sl, exp(thisEta) ) ) * w[w_ind] ;
        con_h = log(lnkFn->con_d(this_d, this_sl, exp(thisEta + this_h) ) ) * w[w_ind] ;
        con_l = log(lnkFn->con_d(this_d, this_sl, exp(thisEta - this_h) ) ) * w[w_ind] ;
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
        
    }
    thisSize = gic.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta = eta[gic[i].nu];
        this_sl = s_v[gic[i].l];
        this_sr = s_v[gic[i].r];
        thisExpEta = exp(thisEta);
        con0 = log(lnkFn->con_s(this_sl, thisExpEta)
                   -lnkFn->con_s(this_sr, thisExpEta) ) * w[w_ind];
        
        thisExpEta = exp(thisEta + this_h);
        
        con_h = log(lnkFn->con_s(this_sl, thisExpEta)
                    -lnkFn->con_s(this_sr, thisExpEta) ) * w[w_ind];
        
        thisExpEta = exp(thisEta - this_h);
        con_l = log(lnkFn->con_s(this_sl, thisExpEta)
                    -lnkFn->con_s(this_sr, thisExpEta) ) * w[w_ind];

        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
        
    }
    thisSize = lc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta = eta[lc[i].nu];
        this_sl = s_v[lc[i].r];
        con0 = log(1.0 - lnkFn->con_s(this_sl, exp(thisEta) ) ) * w[w_ind];        
        con_h = log(1.0 - lnkFn->con_s(this_sl, exp(thisEta + this_h) ) ) * w[w_ind];
        con_l = log(1.0 - lnkFn->con_s(this_sl, exp(thisEta - this_h) ) ) * w[w_ind];
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);

    }
    thisSize = rc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;        
        thisEta = eta[rc[i].nu];
        this_sr = s_v[rc[i].l];
        con0 = log(lnkFn->con_s(this_sr, exp(thisEta) ) ) * w[w_ind];
        con_h = log(lnkFn->con_s(this_sr, exp(thisEta + this_h) ) ) * w[w_ind];
        con_l = log(lnkFn->con_s(this_sr, exp(thisEta - this_h) ) ) * w[w_ind];
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
    }
  
}




void parBLInfo::update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                                     Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals,
                                     Eigen::VectorXd &par){
    for(int i = 0; i < s_t.size(); i++){s_vals[i] = base_s(s_t[i], par);}
    for(int i = 0; i < d_t.size(); i++){d_vals[i] = base_d(d_t[i], par);}
}

void IC_parOpt::calc_baseline_dervs(){
    int k = b_pars.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    d_b_pars.resize(k);
    d2_b_pars.resize(k,k);

    double lk_0 = calcLike_all();
    double org_h = h;
    bool bad_derv = true;
    int tries = 0;
    while(tries < 4 && bad_derv){
        bad_derv = false;
        tries++;
        for(int i = 0; i < k; i++){
            b_pars[i] += h;
            lk_h[i] = calcLike_all();
            b_pars[i] -= 2.0 * h;
            lk_l[i] = calcLike_all();
            b_pars[i] += h;
            d_b_pars[i] = (lk_h[i] - lk_l[i])/(2.0 * h);
            d2_b_pars(i,i) = (lk_h[i] + lk_l[i] - 2.0*lk_0) / (h*h);
        
            if(lk_h[i] == R_NegInf || lk_l[i] == R_NegInf){
                bad_derv = true;
                h = h/4;
            }
        }
    }

    if(bad_derv){Rprintf("error: was not able to calculate derivative of baseline parameters!\n");}
        
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                b_pars[i] += h;
                b_pars[j] += h;
                lk_hh = calcLike_all();
                b_pars[i] -= 2.0 * h;
                b_pars[j] -= 2.0 * h;
                lk_ll = calcLike_all();
                b_pars[i] += h;
                b_pars[j] += h;
                rho = (lk_hh + lk_ll + 2.0 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2.0 * h * h);
                d2_b_pars(i,j) = rho;
                d2_b_pars(j,i) = rho;
            }
        }
    }
    
    calculate_baseline_probs();
    h = org_h;
}


void IC_parOpt::NR_baseline_pars(){
    calc_baseline_dervs();

    double lk_0 = calcLike_baseReady();
    int k = b_pars.size();
    Eigen::VectorXd propVec(k);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve(d2_b_pars);
    Eigen::VectorXd evals(1);
    evals[0] = 1;

    if(esolve.info() == Eigen::Success)
        evals = esolve.eigenvalues();
    
    if(max(evals) < -0.001) { propVec = -d2_b_pars.ldlt().solve(d_b_pars); }
    else{
        for(int i = 0; i < k; i++){
                if(d2_b_pars(i,i) < -0.001) propVec[i] = -d_b_pars[i]/d2_b_pars(i,i);
                else propVec[i] = signVal(d_b_pars[i]) * 0.001;
            
                if(ISNAN(propVec[i])) propVec[i] = 0;
        }
    }
    
    int tries = 0;
    b_pars += propVec;
    propVec *= -1;
    double lk_new = calcLike_all();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        b_pars += propVec;
        lk_new = calcLike_all();
    }
    if(lk_new < lk_0){
        b_pars += propVec;
        lk_new = calcLike_all();
    }
}


void IC_parOpt::update_etas(){
    eta = covars * betas;
    for(int i = 0; i < eta.size(); i++)
        expEta[i] = exp(eta[i]);
}

void IC_parOpt::partAnalyticCovar_dervs(){
    update_dobs_detas();
    
    int n = eta.size();
    int k = betas.size();
    d_betas.resize(k);
    d2_betas.resize(k,k);

    for(int i = 0; i < k; i++){
        d_betas[i] = 0;
        for(int j = 0; j < k; j++){
            d2_betas(i, j) = 0;
        }
    }
    
    double this_d1;
    double this_d2;
    double this_cov_ij;
    for(int i = 0; i < n; i++){
        this_d1 = dobs_deta[i];
        this_d2 = d2obs_d2eta[i];
        for(int j = 0; j < k; j++){
            this_cov_ij = covars(i,j);
            d_betas[j] += this_d1 * this_cov_ij;
            for(int m = 0; m <=j ; m++){
                d2_betas(m ,j) += this_d2 * this_cov_ij * covars(i, m) ;
            }
        }
    }

    for(int i = 0; i < k; i++){
        for(int j = i + 1; j < k; j++){ d2_betas(j,i) = d2_betas(i, j); }
    }
}

void IC_parOpt::numericCovar_dervs(){
    int k = betas.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    d_betas.resize(k);
    d2_betas.resize(k,k);
    
    double lk_0 = calcLike_baseReady();
    
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            d2_betas(i,j) = 0;
        }
    }
    
    
    for(int i = 0; i < k; i++){
        betas[i] += h;
        update_etas();
        lk_h[i] = calcLike_baseReady();
        betas[i] -= 2.0 * h;
        update_etas();
        lk_l[i] = calcLike_baseReady();
        betas[i] += h;
        d_betas[i] = (lk_h[i] - lk_l[i])/(2.0 * h);
        d2_betas(i,i) = (lk_h[i] + lk_l[i] - 2.0 *lk_0) / (h*h);
    }
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                betas[i] += h;
                betas[j] += h;
                update_etas();
                lk_hh = calcLike_baseReady();
                betas[i] -= 2.0 * h;
                betas[j] -= 2.0 * h;
                update_etas();
                lk_ll = calcLike_baseReady();
                betas[i] += h;
                betas[j] += h;
                rho = (lk_hh + lk_ll + 2.0 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2.0 * h * h);
                d2_betas(i,j) = rho;
                d2_betas(j,i) = rho;
            }
        }
    }
    update_etas();
}


void IC_parOpt::fillFullHessianAndScore(SEXP r_mat, SEXP score){
    int k_base = b_pars.size();
    int k_reg = betas.size();
    int k_tot = k_base + k_reg;
    double lk_0 = calcLike_all();
    vector<double> lk_l(k_tot);
    vector<double> lk_h(k_tot);
    for(int i = 0; i < k_base; i++){
        b_pars[i] += h;
        lk_h[i] = calcLike_all();
        b_pars[i] -= 2.0 * h;
        lk_l[i] = calcLike_all();
        b_pars[i] += h;
        REAL(r_mat)[i + i * k_tot] = (lk_h[i] + lk_l[i] - 2.0 * lk_0)/(h*h);
        REAL(score)[i] = (lk_h[i] - lk_l[i])/(2.0*h);
    }
    calculate_baseline_probs();
    int i_tot;
    for(int i = 0; i < k_reg; i++){
        i_tot = i + k_base;
        betas[i] += h;
        update_etas();
        lk_h[i_tot] = calcLike_baseReady();
        betas[i] -= 2.0 * h;
        update_etas();
        lk_l[i_tot] = calcLike_baseReady();
        betas[i] += h;
        REAL(r_mat)[i_tot + i_tot * k_tot] = (lk_l[i_tot] + lk_h[i_tot] - 2.0 * lk_0)/(h * h);
        REAL(score)[i_tot] = (lk_h[i_tot] - lk_l[i_tot])/(2.0 * h);

    }
    update_etas();
    partAnalyticCovar_dervs();
    
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k_tot; i++){
        for(int j = 0; j < i; j++){
            if(i < k_base || j < k_base){
                if(i <k_base)   {   b_pars[i] += h;}
                else            {   betas[i - k_base] += h;}
                
                if(j < k_base)  {   b_pars[j] += h;}
                else            {   betas[j - k_base] += h;};
                update_etas();
                lk_hh = calcLike_all();

                if(i <k_base)   {   b_pars[i] -= 2.0 * h;}
                else            {   betas[i - k_base] -= 2.0 * h;}
                
                if(j < k_base)  {   b_pars[j] -= 2.0 * h;}
                else            {   betas[j - k_base] -= 2.0 * h;};
                update_etas();
                lk_ll = calcLike_all();

                if(i <k_base)   {   b_pars[i] += h;}
                else            {   betas[i - k_base] += h;}
                
                if(j < k_base)  {   b_pars[j] += h;}
                else            {   betas[j - k_base] += h;};
                
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                REAL(r_mat)[i + j * k_tot] = rho;
                REAL(r_mat)[j + i * k_tot] = rho;
            }
            else{
                REAL(r_mat)[i + j * k_tot] = d2_betas(i - k_base, j - k_base);
                REAL(r_mat)[j + i * k_tot] = d2_betas(i - k_base, j - k_base);
            }
        }
    }
    update_etas();
    calculate_baseline_probs();

}

void IC_parOpt::NR_reg_pars(){
    int k = betas.size();
    if(k == 0) return;
    
    partAnalyticCovar_dervs();
//    numericCovar_dervs();
    
    double lk_0 = calcLike_baseReady();

    Eigen::VectorXd propVec(k);
    propVec = -d2_betas.fullPivLu().solve(d_betas);
    
    double err = (d2_betas*propVec + d_betas).norm() / d_betas.norm();
    if(err > .001){
        for(int i = 0; i < k; i++){
            propVec[i] = 0;
            if(d2_betas(i,i) < 0)   propVec[i] = -d_betas[i] / d2_betas(i,i);
            else propVec[i] = signVal(d_betas[i]) * 0.01;
            if(ISNAN(propVec[i])) propVec[i] = 0;
        }
    }
    int tries = 0;
    betas += propVec;
    propVec *= -1;
    update_etas();
    double lk_new = calcLike_all();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        betas += propVec;
        update_etas();
        lk_new = calcLike_all();
    }
    if(lk_new < lk_0){
        betas += propVec;
        update_etas();
        lk_new = calcLike_all();
    }
}


IC_parOpt::IC_parOpt(Rcpp::List R_list){

	Rcpp::NumericVector R_s_t 	   = R_list["s_t"];
	Rcpp::NumericVector R_d_t      = R_list["d_t"];
	Rcpp::NumericMatrix R_covars   = R_list["covars"];
	Rcpp::IntegerMatrix R_uncenInd = R_list["uncenInd_mat"];
	Rcpp::IntegerMatrix R_gicInd   = R_list["gicInd_mat"];
	Rcpp::IntegerVector R_lInd     = R_list["leftCenInd"];
	Rcpp::IntegerVector R_rInd     = R_list["rightCenInd"];
	Rcpp::IntegerVector R_parType  = R_list["parInd"];
	Rcpp::IntegerVector R_linkType = R_list["linkType"];
	Rcpp::NumericVector R_w        = R_list["w"];
	
	successfulBuild = true;
	if(Rf_isNull(R_s_t) )     successfulBuild = false;
	if(Rf_isNull(R_d_t) )     successfulBuild = false;
	if(Rf_isNull(R_covars))   successfulBuild = false;
	if(Rf_isNull(R_uncenInd)) successfulBuild = false;
	if(Rf_isNull(R_gicInd))   successfulBuild = false;
	if(Rf_isNull(R_lInd))     successfulBuild = false;
	if(Rf_isNull(R_parType))  successfulBuild = false;
	if(Rf_isNull(R_linkType)) successfulBuild = false;
	if(Rf_isNull(R_w))        successfulBuild = false;
	
	
	if(!successfulBuild){
		Rprintf("Build unsuccessful because list names are not correct!\n");
	}
	
	init(R_s_t, R_d_t, R_covars, 
 		 R_uncenInd, R_gicInd, R_lInd,
 		 R_rInd, R_parType, R_linkType, 
 		 R_w);
}

IC_parOpt::IC_parOpt(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
                     SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
                     SEXP R_parType, SEXP R_linkType, SEXP R_w){
 	init(R_s_t, R_d_t, R_covars, 
 		 R_uncenInd, R_gicInd, R_lInd,
 		 R_rInd, R_parType, R_linkType, 
 		 R_w);                    
}

void IC_parOpt::init(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
                     SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
                     SEXP R_parType, SEXP R_linkType, SEXP R_w){
    blInf = NULL;
    parType = INTEGER(R_parType)[0];
    if(INTEGER(R_parType)[0] == 1) {
        blInf = new gammaInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 2){
        blInf = new weibullInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 3){
        blInf = new lnormInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 4){
        blInf = new expInfo();
        b_pars.resize(1);
        b_pars[0] = 0;
    }
    else if(INTEGER(R_parType)[0] == 5){
        blInf = new loglogisticInfo();
        b_pars.resize(2);
        b_pars[0] = 0;
        b_pars[1] = 0;
    }
    else if(INTEGER(R_parType)[0] == 6){
        blInf = new genGammaInfo();
        b_pars.resize(3);
        b_pars[0] = 0;
        b_pars[1] = 0;
        b_pars[2] = 0;
    }
    else{Rprintf("warning: parameter type not supported!\n");}
    
    lnkFn = NULL;
    linkType = INTEGER(R_linkType)[0];
    if(INTEGER(R_linkType)[0] == 1) {lnkFn = new propOdd;}
    else if(INTEGER(R_linkType)[0] == 2) {lnkFn = new propHaz;}
    else if(INTEGER(R_linkType)[0] == 3) {lnkFn = new aft_linkFun;}
    else{Rprintf("warning: link type not supported!\n");}
    
    Rvec2eigen(R_s_t, s_t);
    Rvec2eigen(R_d_t, d_t);
    s_v.resize(s_t.size());
    d_v.resize(d_t.size());
    copyRmatrix_intoEigen(R_covars, covars);
    int k = covars.cols();
    betas.resize(k);
    for(int i = 0; i < k; i++)  betas[i] = 0;
    d_betas.resize(k);
    d2_betas.resize(k, k);
    
    SEXP RuncenDim = Rf_getAttrib(R_uncenInd, R_DimSymbol);
    PROTECT(RuncenDim);
    SEXP RgicDim = Rf_getAttrib(R_gicInd, R_DimSymbol);
    PROTECT(RgicDim);
    
    int n_1 = INTEGER(RuncenDim)[0];
    int n_2 = INTEGER(RgicDim)[0];
    int n_3 = LENGTH(R_lInd);
    int n_4 = LENGTH(R_rInd);
    
    int tot_n = n_1 + n_2 + n_3 + n_4;
    eta.resize(tot_n);
    expEta.resize(tot_n);
    w.resize(tot_n);
    dobs_deta.resize(tot_n);
    d2obs_d2eta.resize(tot_n);
    
    for(int i = 0; i < tot_n; i++){
        eta[i] = 0;
        expEta[i] = 1;
        w[i] = REAL(R_w)[i];
        dobs_deta[i] = 0;
        d2obs_d2eta[i] = 0;
    }
    
    uc.resize(n_1);
    for(int i = 0; i < n_1; i++){
        uc[i].d = INTEGER(R_uncenInd)[i] - 1;
        uc[i].s = INTEGER(R_uncenInd)[i + n_1] - 1;
        uc[i].nu = i;
    }
    
    gic.resize(n_2);
    for(int i = 0; i < n_2; i++){
        gic[i].l = INTEGER(R_gicInd)[i] - 1;
        gic[i].r = INTEGER(R_gicInd)[i + n_2] - 1;
        gic[i].nu = i + n_1;
    }
    
    lc.resize(n_3);
    for(int i = 0; i < n_3; i++){
        lc[i].r = INTEGER(R_lInd)[i] - 1;
        lc[i].nu = i + n_1 + n_2;
    }
    
    rc.resize(n_4);
    for(int i = 0; i < n_4; i++){
        rc[i].l = INTEGER(R_rInd)[i] - 1;
        rc[i].nu = i + n_1 + n_2 + n_3;
    }
    
    h = pow(10.0, -5.0);
    UNPROTECT(2);
    iter = 0;
    
    successfulBuild = true;
}

void IC_parOpt::optimize(){
    lk_old = R_NegInf;
    int maxIter = 1000;
    double tol = pow(10.0, -10.0);
    lk_new = calcLike_all();

    if(lk_new == R_NegInf){
        int bk = b_pars.size();
        int tries = 0;
        double delta = 0.001;
        while(tries < 10 && lk_new == R_NegInf){
            tries++;
            for(int i = 0; i < bk; i++){
                if(lk_new == R_NegInf){
                    b_pars[i] = delta;
                    lk_new = calcLike_all();
                    if(lk_new == R_NegInf) b_pars[i] = 0;
                }
            }
            delta *= 5.0;
        }
    }
    if(lk_new == R_NegInf){
        int bk = b_pars.size();
        int tries = 0;
        double delta = -1.;
        while(tries < 10 && lk_new == R_NegInf){
            tries++;
            for(int i = 0; i < bk; i++){
                if(lk_new == R_NegInf){
                    b_pars[i] = delta;
                    lk_new = calcLike_all();
                    if(lk_new == R_NegInf) b_pars[i] = 0;
                }
            }
            delta *= 5.0;
        }
    }
    if(lk_new == R_NegInf){
        Rprintf("failed to find adequate starting point!n");
        return;
    }

    for(int i = 0; i < 5; i++){ NR_baseline_pars(); }
    while(iter < maxIter && lk_new - lk_old > tol){
        lk_old = lk_new;
        iter++;
        NR_baseline_pars();
        NR_reg_pars();
        lk_new = calcLike_baseReady();   
    }
        
}

Rcpp::List ic_par(SEXP R_s_t, SEXP R_d_t, SEXP covars,
            SEXP uncenInd, SEXP gicInd, SEXP lInd, SEXP rInd,
            SEXP parType, SEXP linkType,
            SEXP outHessian, SEXP R_w){
    IC_parOpt* optObj;
    if(INTEGER(linkType)[0] == 1 || INTEGER(linkType)[0] == 2){
    	optObj = new IC_parOpt(R_s_t, R_d_t, covars, uncenInd, gicInd, 
    						   lInd, rInd, parType, linkType, R_w);
    }
    else if(INTEGER(linkType)[0] == 3){
    	optObj = new IC_parOpt_aft(R_s_t, R_d_t, covars, uncenInd, gicInd, 
    						   lInd, rInd, parType, linkType, R_w);
    	
    }
    else{
    	Rprintf("Warning: linkType not recognized.\n");
    	return(R_NilValue);
    }
    if(optObj->blInf == NULL) return(R_NilValue);
    if(optObj->lnkFn == NULL) return(R_NilValue);

	optObj->optimize();
	Rcpp::List ans;
	ans = optObj->exportAns();
	delete optObj;
	return(ans);
}


Rcpp::List ic_parList(Rcpp::List R_parList){
    IC_parOpt* optObj;
    
    Rcpp::IntegerVector linkType = R_parList["linkType"];
        
    if(INTEGER(linkType)[0] == 1 || INTEGER(linkType)[0] == 2){
    	optObj = new IC_parOpt(R_parList);
    }
    else if(INTEGER(linkType)[0] == 3){
    	optObj = new IC_parOpt_aft(R_parList);
    	
    }
    else{
    	Rprintf("Warning: linkType not recognized.\n");
    	return(R_NilValue);
    }


    if(optObj->blInf == NULL) return(R_NilValue);
    if(optObj->lnkFn == NULL) return(R_NilValue);

	optObj->optimize();
	Rcpp::List ans = optObj->exportAns();
	delete optObj;
	return(ans);
}

Rcpp::List IC_parOpt::exportAns(){
	int totParams = betas.size() + b_pars.size();

    Rcpp::NumericMatrix outHessian(totParams, totParams);
	Rcpp::NumericVector score(totParams);
	Rcpp::NumericVector reg_est( betas.size() );
	Rcpp::NumericVector base_est( b_pars.size() );
	Rcpp::NumericVector final_llk(1);	
	Rcpp::NumericVector iters(1);
		
	fillFullHessianAndScore(outHessian, score);
    for(int i = 0; i < LENGTH(reg_est); i++)    reg_est[i] = betas[i];
    for(int i = 0; i < LENGTH(base_est); i++)   base_est[i] = b_pars[i];
    final_llk[0] = calcLike_baseReady();
    iters[0]     = iter;
        
    Rcpp::List ans(6);
    ans["reg_pars"]   = reg_est;
    ans["baseline"]   = base_est;
    ans["llk"]        = final_llk;
    ans["iterations"] = iters;
    ans["hessian"]    = outHessian;
    ans["score"]      = score;

	return(ans);
}

IC_parOpt::~IC_parOpt(){
    if(parType == 1){
        gammaInfo* deleteObj = static_cast<gammaInfo*>(blInf);
        delete deleteObj;
    }
    if(parType == 2){
        weibullInfo* deleteObj = static_cast<weibullInfo*>(blInf);
        delete deleteObj;
    }
    if(parType == 3){
        lnormInfo* deleteObj = static_cast<lnormInfo*>(blInf);
        delete deleteObj;
    }
    if(parType == 4){
        expInfo* deleteObj = static_cast<expInfo*>(blInf);
        delete deleteObj;
    }
    if(parType == 5){
        loglogisticInfo* deleteObj = static_cast<loglogisticInfo*>(blInf);
        delete deleteObj;
    }
    if(parType == 6){
        genGammaInfo* deleteObj = static_cast<genGammaInfo*>(blInf);
        delete deleteObj;
    }
    if(linkType == 1){
        propOdd* deleteObj = static_cast<propOdd*>(lnkFn);
        delete deleteObj;
//	    delete optObj;
    }
    if(linkType == 2){
        propHaz* deleteObj = static_cast<propHaz*>(lnkFn);
        delete deleteObj;
//	    delete optObj;
    }
    if(linkType == 3){
    	aft_linkFun* deleteLnkObj = static_cast<aft_linkFun*>(lnkFn);
    	delete deleteLnkObj;
//    	IC_parOpt_aft* deleteObj = static_cast<IC_parOpt_aft*>(optObj);
//    	delete deleteObj;
    }
}








SEXP dGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q){
    int size = LENGTH(R_x);
    
    double* x = REAL(R_x);
    double* mu = REAL(R_mu);
    double* s = REAL(R_s);
    double* Q = REAL(R_Q);
    
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, size));
    double* cans = REAL(ans);
    for(int i = 0; i < size; i++){
        cans[i] = ic_dgeneralgamma(x[i], mu[i], s[i], Q[i]);
    }
    UNPROTECT(1);
    return(ans);
}


SEXP pGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q){
    int size = LENGTH(R_x);
    
    double* x = REAL(R_x);
    double* mu = REAL(R_mu);
    double* s = REAL(R_s);
    double* Q = REAL(R_Q);
    
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, size));
    double* cans = REAL(ans);
    for(int i = 0; i < size; i++){
        cans[i] = ic_pgeneralgamma(x[i], mu[i], s[i], Q[i]);
    }
    UNPROTECT(1);
    return(ans);
}
SEXP qGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q){
    int size = LENGTH(R_x);
    double* x = REAL(R_x);
    double* mu = REAL(R_mu);
    double* s = REAL(R_s);
    double* Q = REAL(R_Q);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP, size));
    double* cans = REAL(ans);
    for(int i = 0; i < size; i++){
        cans[i] = ic_qgeneralgamma(x[i], mu[i], s[i], Q[i]);
    }
    UNPROTECT(1);
    return(ans);
}
