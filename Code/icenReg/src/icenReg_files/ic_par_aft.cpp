
double IC_parOpt_aft::calcLike_baseReady(){
  double ans = 0;
    int w_ind = -1;
    int thisSize = uc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(con_d(d_t[uc[i].d], expEta[uc[i].nu])) * w[w_ind];
    }
    thisSize = gic.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(con_s(s_t[gic[i].l], expEta[gic[i].nu])
                   -con_s(s_t[gic[i].r], expEta[gic[i].nu]) ) * w[w_ind];
    }
    thisSize = lc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(1.0 - con_s(s_t[lc[i].r], expEta[lc[i].nu])) * w[w_ind];
    }
    thisSize = rc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        ans += log(con_s(s_t[rc[i].l], expEta[rc[i].nu])) * w[w_ind];
    }
    
    if(ISNAN(ans)) ans = R_NegInf;
    return(ans);
}

void IC_parOpt_aft::update_dobs_detas(){
    double con0, con_l, con_h, thisEta, thisExpEta, this_d, this_sl, this_sr;
    int w_ind = -1;
    int thisSize = uc.size();
    
    double this_h = h * 0.1;
    
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta  = eta[uc[i].nu];
        this_d = d_t[uc[i].d];
        con0 = log(con_d(this_d, exp(thisEta) ) ) * w[w_ind] ;
        con_h = log(con_d(this_d, exp(thisEta + this_h) ) ) * w[w_ind] ;
        con_l = log(con_d(this_d, exp(thisEta - this_h) ) ) * w[w_ind] ;
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
        
    }
    thisSize = gic.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta = eta[gic[i].nu];
        this_sl = s_t[gic[i].l];
        this_sr = s_t[gic[i].r];
        thisExpEta = exp(thisEta);
        con0 = log(con_s(this_sl, thisExpEta)
                   -con_s(this_sr, thisExpEta) ) * w[w_ind];
        
        thisExpEta = exp(thisEta + this_h);
        
        con_h = log(con_s(this_sl, thisExpEta)
                    -con_s(this_sr, thisExpEta) ) * w[w_ind];
        
        thisExpEta = exp(thisEta - this_h);
        con_l = log(con_s(this_sl, thisExpEta)
                    -con_s(this_sr, thisExpEta) ) * w[w_ind];

        thisExpEta = exp(thisEta + 2.0 * this_h);
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
        
    }
    thisSize = lc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;
        thisEta = eta[lc[i].nu];
        this_sl = s_t[lc[i].r];
        con0 = log(1.0 - con_s(this_sl, exp(thisEta) ) ) * w[w_ind];        
        con_h = log(1.0 - con_s(this_sl, exp(thisEta + this_h) ) ) * w[w_ind];
        con_l = log(1.0 - con_s(this_sl, exp(thisEta - this_h) ) ) * w[w_ind];
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);

    }
    thisSize = rc.size();
    for(int i = 0; i < thisSize; i++){
        w_ind++;        
        thisEta = eta[rc[i].nu];
        this_sr = s_t[rc[i].l];
        con0 = log(con_s(this_sr, exp(thisEta) ) ) * w[w_ind];
        con_h = log(con_s(this_sr, exp(thisEta + this_h) ) ) * w[w_ind];
        con_l = log(con_s(this_sr, exp(thisEta - this_h) ) ) * w[w_ind];
        dobs_deta[w_ind] = (con_h - con_l) / (2.0 * this_h);
        d2obs_d2eta[w_ind] = (con_h + con_l - 2.0 * con0) / (this_h * this_h);
    }
  
}

IC_parOpt_aft::IC_parOpt_aft(Rcpp::List R_list){
	
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
	
	init(R_s_t, R_d_t, R_covars, 
 		 R_uncenInd, R_gicInd, R_lInd,
 		 R_rInd, R_parType, R_linkType, 
 		 R_w);

}

IC_parOpt_aft::IC_parOpt_aft(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
                     SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
                     SEXP R_parType, SEXP R_linkType, SEXP R_w){
    blInf = NULL;
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
}
