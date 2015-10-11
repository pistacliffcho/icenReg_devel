//
//  ic_sp_cm.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#include "ic_sp_ch.h"
//#include "../icenReg_files/basicUtilities.cpp"


/*      LIKELIHOOD TOOLS        */
void icm_Abst::update_p_ob(int i){
    double chl = baseCH[ obs_inf[i].l ];
    double chr = baseCH[ obs_inf[i].r +1 ];
    double eta = etas[i];
    obs_inf[i].pob = basHaz2CondS(chl, eta) - basHaz2CondS(chr, eta);
}

double icm_Abst::sum_llk(){
    int n = obs_inf.size();
    double ans = 0;
    for(int i = 0; i < n; i++){
        update_p_ob(i);
        ans += log(obs_inf[i].pob) * w[i];
    }
    if(ISNAN(ans)) {ans = R_NegInf;}
    return(ans);
}

double icm_Abst::par_llk(int ind){
    int num_l = node_inf[ind].l.size();
    int num_r = node_inf[ind].r.size();
    double ans = 0;
    int thisInd;
    for(int i = 0; i < num_l; i++){
        thisInd = node_inf[ind].l[i];
        update_p_ob(thisInd);
        ans+= log(obs_inf[thisInd].pob) * w[thisInd];
    }
    for(int i = 0; i < num_r; i++){
        thisInd = node_inf[ind].r[i];
        update_p_ob(thisInd);
        ans+= log(obs_inf[thisInd].pob) * w[thisInd];
    }
    if(ISNAN(ans)) ans = R_NegInf;
    return(ans);
}

void icm_Abst::update_etas(){
    etas = covars * reg_par;
    for(int i = 0; i < etas.size(); i++)
        expEtas[i] = exp(etas[i]);
}

void cumhaz2p_hat(Eigen::VectorXd &ch, vector<double> &p){
    int k = ch.size();
    vector<double> S(k);
    p.resize(k-1);
    for(int i = 0; i < k; i++){ S[i] = exp(-exp(ch[i])); }
    for(int i = 0; i < (k-1); i++){ p[i] = S[i+1] - S[i]; }
}


void icm_Abst::icm_addPar(vector<double> &delta){
    int p_k = delta.size();
    int a_k = baseCH.size();
    if( (p_k+2) != a_k){Rprintf("in icm_addPar, delta is not the same length as actIndex!\n");return;}
    for(int i = 0; i < p_k; i++){ baseCH[i+1] += delta[i]; }
}




/*      INITIALIZATION TOOLS    */
void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, icm_Abst* icm_obj){
    icm_obj->h = 0.0001;
    int n = LENGTH(Rlind);
    if(n != LENGTH(Rrind)){Rprintf("length of Rlind and Rrind not equal\n"); return;}
    icm_obj->base_p_obs.resize(n);
    icm_obj->etas.resize(n);
    icm_obj->expEtas.resize(n);
    icm_obj->w.resize(n);
    
    for(int i = 0; i < n; i++){
        icm_obj->etas[i]       = 0;
        icm_obj->expEtas[i]    = 1;
        icm_obj->base_p_obs[i] = 0;
        icm_obj->w[i]          = REAL(R_w)[i];
    }
    
    copyRmatrix_intoEigen(RCovars, icm_obj->covars);
    int reg_k = icm_obj->covars.cols();
    if(reg_k == 0) icm_obj->hasCovars = false; else icm_obj->hasCovars = true;
    if(reg_k > 0){
        if(n != icm_obj->covars.rows()) {Rprintf("covar rows not equal to n!\n"); return;}
    }
    icm_obj->reg_d1.resize(reg_k);
//    icm_obj->reg_d2.resize(reg_k, reg_k);
    icm_obj->reg_d2.resize(reg_k);
    icm_obj->reg_par.resize(reg_k);
    for(int i = 0; i < reg_k; i++){ icm_obj->reg_par[i] = 0; }
    
    
    int maxInd = 0;
    for(int i = 0; i < n; i++){
        maxInd = max(maxInd, INTEGER(Rrind)[i]);
    }

    

    icm_obj->baseCH.resize(maxInd + 2);
    for(int i = 0; i <= maxInd; i++){ icm_obj->baseCH[i] = R_NegInf; }
    icm_obj->baseCH[maxInd+1] = R_PosInf;
    
/*    icm_obj->H_d1.resize(maxInd - 1);
    icm_obj->H_d2.resize(maxInd - 1, maxInd - 1);   */
    
    vector<int> minActPoints;
    int this_l, this_r;
    icm_obj->obs_inf.resize(n);
    icm_obj->node_inf.resize(maxInd + 2);
    
    
    for(int i = 0; i < n; i++){
        this_l = INTEGER(Rlind)[i];
        this_r = INTEGER(Rrind)[i];
        addIfNeeded(minActPoints, this_l, this_r, maxInd);
        icm_obj->obs_inf[i].l = this_l;
        icm_obj->obs_inf[i].r = this_r;
        icm_obj->node_inf[this_l].l.push_back(i);
        icm_obj->node_inf[this_r + 1].r.push_back(i);
    }
    
    sort(minActPoints.begin(), minActPoints.end());
    
    
    double stepSize = 4.0/minActPoints.size();
    double curVal = -2.0;
    
    int intActIndex = 0;
    for(int i = 1; i < (maxInd+1); i++){
        if(intActIndex < minActPoints.size() ){
            if(i == minActPoints[intActIndex]){
                intActIndex++;
                curVal += stepSize;
            }
        }
        icm_obj->baseCH[i] = curVal;
    }
    icm_obj->startGD = false;
//    double gd_mult = 1.0;
}



/*      OPTIMIZATION TOOLS      */
void icm_Abst::numericBaseDervsOne(int raw_ind, vector<double> &dvec){
    dvec.resize(2);
    dvec[0] = 0;
    dvec[1] = 0;
    if(raw_ind <= 0 || raw_ind >= (baseCH.size()- 1)){Rprintf("warning: inappropriate choice of ind for numericBaseDervs ind = %d\n", raw_ind); return;}
    
    baseCH[raw_ind] += h;
    double llk_h = par_llk(raw_ind);
    baseCH[raw_ind] -= 2*h;
    double llk_l = par_llk(raw_ind);
    baseCH[raw_ind] += h;
    double llk_st = par_llk(raw_ind);
    
    dvec[0] = (llk_h - llk_l)/(2*h);
    dvec[1] = (llk_h + llk_l - 2 * llk_st) / (h * h);
}

void icm_Abst::numericBaseDervsAllAct(vector<double> &d1, vector<double> &d2){
    int k = baseCH.size();
    d1.resize(k);
    d2.resize(k);
    vector<double> ind_dervs(2);
    for(int i = 1; i < (k-1); i++){
        numericBaseDervsOne(i, ind_dervs);
        d1[i] = ind_dervs[0];
        d2[i] = ind_dervs[1];
    }
}

void icm_Abst::numericBaseDervsAllRaw(vector<double> &d1, vector<double> &d2){
    int k = baseCH.size() - 2;
    d1.resize(k);
    d2.resize(k);
    vector<double> ind_dervs(2);
    for(int i = 0; i < k; i++){
        numericBaseDervsOne(i + 1, ind_dervs);
        d1[i] = ind_dervs[0];
        d2[i] = ind_dervs[1];
    }
}

 
void icm_Abst::icm_step(){
    
    vector<double> d1;
    vector<double> d2;
    numericBaseDervsAllRaw(d1, d2);
    int thisSize = d1.size();
    for(int i = 0; i < thisSize; i ++){
        if(ISNAN(d2[i]))    {Rprintf("warning: d2 isnan!\n"); return;}
        if(d2[i] >= 0 || d2[i] == R_NegInf)      {Rprintf("warning: invalid d2 in icm step. i = %d, d2 = %f. Quiting icm step\n", i, d2[i]); return;}
    }
    vector<double> x(d1.size());
    if(x.size() != baseCH.size() - 2){Rprintf("warning: x.size()! = actIndex.size()\n"); return;}
    thisSize = baseCH.size() - 2;
    for(int i = 0; i < thisSize; i++){x[i] = baseCH[i + 1];}
    vector<double> prop(d1.size());
    
    double llk_st = sum_llk();
    
    pavaForOptim(d1, d2, x, prop);
    
    icm_addPar(prop);
    double llk_new = sum_llk();
    mult_vec(-1.0, prop);
    int tries = 0;
    while(llk_st > llk_new && tries < 5){
        tries++;
        mult_vec(0.5, prop);
        icm_addPar(prop);
        llk_new = sum_llk();
    }
    if(llk_new < llk_st){
        icm_addPar(prop);
        llk_new = sum_llk();
        mult_vec(0, prop);
    }
    
    maxBaseChg = 0;
    for(int i = 0; i < thisSize; i++){
        maxBaseChg = max(maxBaseChg, abs(prop[i]) );
    }

}

void icm_Abst::calcAnalyticRegDervs(Eigen::VectorXd &hess, Eigen::VectorXd &d1){
    int k = reg_par.size();
    int n = etas.size();
    
    Eigen::VectorXd l_cont(n);
    Eigen::VectorXd r_cont(n);
    Eigen::VectorXd totCont(n);

    Eigen::VectorXd l_cont2(n);
    Eigen::VectorXd r_cont2(n);
    Eigen::VectorXd totCont2(n);
    
    int lind, rind;
    double l_ch, r_ch, eta, pob, log_p;
    for(int i = 0; i < n; i++){
        l_cont[i]  = 0;
        r_cont[i]  = 0;
        l_cont2[i] = 0;
        r_cont2[i] = 0;

        lind = obs_inf[i].l;
        rind = obs_inf[i].r;
        pob  = obs_inf[i].pob;
        log_p = log(pob);
        l_ch = baseCH[lind];
        r_ch = baseCH[rind + 1];
        eta  = etas[i];
        if(l_ch > R_NegInf){
            l_cont[i]  = reg_d1_lnk(l_ch, eta, log_p);
            l_cont2[i] = reg_d2_lnk(l_ch, eta, log_p);
        }
        if(r_ch < R_PosInf){
            r_cont[i]  = -reg_d1_lnk(r_ch, eta, log_p);
            r_cont2[i] = -reg_d2_lnk(r_ch, eta, log_p);
        }
        totCont[i] = l_cont[i] + r_cont[i];
        totCont2[i] = l_cont2[i] + r_cont2[i] - totCont[i] * totCont[i];
    }
    
    hess.resize(k);
    d1.resize(k);
    for(int i = 0; i < k; i++){
        d1[i] = 0;
        hess[i] = 0;
    }

    double this_covar;
    double this_w;
    double this_w_covar;
    for(int i = 0; i < n; i++){
        this_w = w[i];
        for(int a = 0; a < k; a++){
            this_covar = covars(i,a);
            this_w_covar = this_w * this_covar;
            d1[a] += this_w_covar * totCont[i];
    /*        for(int b = 0; b < a; b++){
                hess(a,b) += covars(i,a) * covars(i,b) * totCont2[i];
                hess(b,a) = hess(a,b);
            }       ignorable due to PCA    */
            hess[a] += this_w_covar * this_covar * totCont2[i];
        }
    }

}

/*
void icm_Abst::numericRegDervs(){
    int k = reg_par.size();
    vector<double> lk_l(k);
    vector<double> lk_h(k);
    reg_d1.resize(k);
//    reg_d2.resize(k,k);
    reg_d2.resize(k);
    
    double lk_0 = sum_llk();
    
    
    for(int i = 0; i < k; i++){
        reg_par[i] += h;
        update_etas();
        lk_h[i] = sum_llk();
        reg_par[i] -= 2 * h;
        update_etas();
        lk_l[i] = sum_llk();
        reg_par[i] += h;
        reg_d1[i] = (lk_h[i] - lk_l[i])/(2 * h);
        reg_d2(i,i) = (lk_h[i] + lk_l[i] - 2*lk_0) / (h*h);
    }
    
    double lk_ll, lk_hh, rho;
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            if(i != j){
                reg_par[i] += h;
                reg_par[j] += h;
                update_etas();
                lk_hh = sum_llk();
                reg_par[i] -= 2 * h;
                reg_par[j] -= 2 * h;
                update_etas();
                lk_ll = sum_llk();
                reg_par[i] += h;
                reg_par[j] += h;
                rho = (lk_hh + lk_ll + 2 * lk_0 - lk_h[i] - lk_h[j] - lk_l[i] - lk_l[j])/(2 * h * h);
                reg_d2(i,j) = rho;
                reg_d2(j,i) = rho;
            }
        }
    }
    update_etas();
}
*/
 
 
void icm_Abst::covar_nr_step(){
    int k = reg_par.size();
    calcAnalyticRegDervs(reg_d2, reg_d1);
    double lk_0 = sum_llk();
/*    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esolve(reg_d2);
    Eigen::VectorXd evals(1);
    evals[0] = 1;
    
    if(esolve.info() == Eigen::Success)
        evals = esolve.eigenvalues();
    int tries = 0;
    double delta = 1;
    while(max(evals) > -0.000001 && tries < 10){
        tries++;
        for(int i = 0; i < k; i++)
            reg_d2(i,i) -= delta;
        delta *= 2;
        esolve.compute(reg_d2);
        if(esolve.info() == Eigen::Success)
            evals = esolve.eigenvalues();
    }
    
    if(max(evals) > 0){
        return;
    }       */
    for(int i = 0; i < k; i++){
        if(reg_d2[i] >= -0.0000001 || ISNAN(reg_d2[i])){
            reg_d2[i] = -100.00;
        }
        if(ISNAN(reg_d1[i]) ){reg_d1[i] = 0;}
    }
 
//    Eigen::VectorXd propVec = -reg_d2.ldlt().solve(reg_d1);
    Eigen::VectorXd propVec(k);
    for(int i = 0; i < k; i++){propVec[i] = -reg_d1[i]/reg_d2[i];}
    int tries = 0;
    reg_par += propVec;
    propVec *= -1;
    update_etas();
    double lk_new = sum_llk();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        reg_par += propVec;
        update_etas();
        lk_new = sum_llk();
    }
}


/*      CALLING ALGORITHM FROM R     */
SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType, SEXP R_w){
    icm_Abst* optObj;
    if(INTEGER(fitType)[0] == 1){
        optObj = new icm_ph;
    }
    else if(INTEGER(fitType)[0] == 2){
        optObj = new icm_po;
    }
    else { Rprintf("fit type not supported\n");return(R_NilValue);}

    setup_icm(Rlind, Rrind, Rcovars, R_w, optObj);
    
    double gd_llk_1, gd_llk_2;
    
    double llk_old = R_NegInf;
    double llk_new = optObj->sum_llk();
    int tries = 0;
    
    bool metOnce = false;
    double tol = pow(10, -10);
    while(tries < 500 && (llk_new - llk_old) > tol){
        tries++;
//        if(tries % 20 == 0){R_CheckUserInterrupt();}
        llk_old = llk_new;
        if(optObj->hasCovars)       optObj->covar_nr_step();
        for(int i = 0; i < 10; i++)  {
            optObj->icm_step();
            
            gd_llk_1 = optObj->sum_llk();
            
            optObj->gradientDescent_step();

            gd_llk_2 = optObj->sum_llk();
            Rprintf("change from gradientDescent step = %f\n", gd_llk_2 - gd_llk_1);
 
            
            if(optObj->maxBaseChg < pow(10.0, -12.0)){break;}
        }
        llk_new = optObj->sum_llk();
        
        if(llk_new - llk_old > tol){metOnce = false;}
        if(metOnce == false){
            if(llk_new - llk_old <= tol){
                metOnce = true;
                llk_old = llk_old - 2 * tol;
            }
        }
    }

 
    if((llk_new - llk_old) < -0.001 ){
        Rprintf("warning: likelihood decreased! difference = %f\n", llk_new - llk_old);
    }
    
    vector<double> p_hat;
    cumhaz2p_hat(optObj->baseCH, p_hat);
    
    
    SEXP ans = PROTECT(allocVector(VECSXP, 5));
    SEXP R_pans = PROTECT(allocVector(REALSXP,p_hat.size()));
    SEXP R_coef = PROTECT(allocVector(REALSXP, optObj->reg_par.size()));
    SEXP R_fnl_llk = PROTECT(allocVector(REALSXP, 1));
    SEXP R_its = PROTECT(allocVector(REALSXP, 1));
    SEXP R_score = PROTECT(allocVector(REALSXP, optObj->reg_par.size()));
    int phat_size = p_hat.size();
    for(int i = 0; i < phat_size; i++){ REAL(R_pans)[i] = p_hat[i]; }
    for(int i = 0; i < optObj->reg_par.size(); i++){
        REAL(R_coef)[i] = optObj->reg_par[i];
        REAL(R_score)[i] = optObj->reg_d1[i];
    }
    REAL(R_fnl_llk)[0] = llk_new;
    REAL(R_its)[0] = tries;
    
    SET_VECTOR_ELT(ans, 0, R_pans);
    SET_VECTOR_ELT(ans, 1, R_coef);
    SET_VECTOR_ELT(ans, 2, R_fnl_llk);
    SET_VECTOR_ELT(ans, 3, R_its);
    SET_VECTOR_ELT(ans, 4, R_score);
    
    UNPROTECT(6);
    
    if(INTEGER(fitType)[0] == 1){
        icm_ph* deleteObj = static_cast<icm_ph*>(optObj);
        delete deleteObj;
    }
    else if(INTEGER(fitType)[0] == 2){
        icm_po* deleteObj = static_cast<icm_po*>(optObj);
        delete deleteObj;
    }
    
    return(ans);

}


/*      GETTING MAXIMAL INTERSECTIONS       */

SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals){
    //NOTE: R_AllVals MUST be sorted!!
    int k = LENGTH(R_AllVals);
    vector<double> mi_l;
    vector<double> mi_r;
    
    bool foundLeft = false;
    double last_left = R_NegInf;
    for(int i = 0; i < k; i++){
        if(!foundLeft)                      foundLeft = LOGICAL(isL)[i] == TRUE;
        if(LOGICAL(isL)[i] == TRUE)         last_left = REAL(R_AllVals)[i];
        if(foundLeft){
            if(LOGICAL(isR)[i] == TRUE){
                mi_l.push_back(last_left);
                mi_r.push_back(REAL(R_AllVals)[i]);
                foundLeft = false;
            }
        }
    }
    int tbulls = mi_l.size();
    
    int n = LENGTH(lVals);
    SEXP l_ind = PROTECT(allocVector(INTSXP, n));
    SEXP r_ind = PROTECT(allocVector(INTSXP, n));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < tbulls; j++){
            if(mi_l[j] >= REAL(lVals)[i]){
                INTEGER(l_ind)[i] = j;
                break;
            }
        }
        
        for(int j = tbulls-1; j >= 0; j--){
            if(mi_r[j] <= REAL(rVals)[i]){
                INTEGER(r_ind)[i] = j;
                break;
            }
        }
    }
    
    
    SEXP ans = PROTECT(allocVector(VECSXP, 4));
    SEXP Rl_mi = PROTECT(allocVector(REALSXP, tbulls));
    SEXP Rr_mi = PROTECT(allocVector(REALSXP, tbulls));
    for(int i = 0; i < tbulls; i++){
        REAL(Rl_mi)[i] = mi_l[i];
        REAL(Rr_mi)[i] = mi_r[i];
    }
    SET_VECTOR_ELT(ans, 0, l_ind);
    SET_VECTOR_ELT(ans, 1, r_ind);
    SET_VECTOR_ELT(ans, 2, Rl_mi);
    SET_VECTOR_ELT(ans, 3, Rr_mi);
    UNPROTECT(5);
    return(ans);
}
