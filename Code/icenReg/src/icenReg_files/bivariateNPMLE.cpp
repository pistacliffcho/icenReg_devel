//
//  bivariateNPMLE.cpp
//  
//
//  Created by Cliff Anderson Bergman on 8/29/15.
//
//

#include "bivariateNPMLE.h"

void bvcen::squeeze(int i, int j){
    if(p_mass[i] == 0 && p_mass[j] == 0){return;}
    std::vector<int> in1not2, in2not1;
    findIndexDiffs(pmass_in_ob[i], pmass_in_ob[j], in1not2, in2not1);
    squeezeInternal(i, j, in1not2, in2not1);
}

void bvcen::squeezeInternal(int p1_ind, int p2_ind, std::vector<int> &in1not2, std::vector<int> &in2not1){
   
    
    double p1 = p_mass[p1_ind];
    double p2 = p_mass[p2_ind];

    double new_p1, new_p2;
    
    double b0 = p1 + p2;
    if(b0 == 0){return;}
    double b1 = 2.0;
    double b2 = 2.0;
    
    int k1 = in1not2.size();
    int k2 = in2not1.size();
    
    int thisInd;

    
/*    if(k1 == 0){
        double new_p1 = 0;
        double new_p2 = p1 + p2;
        
        double d_p1 = -p1;
        double d_p2 = p1;
        for(int i = 0; i < k1; i++){
            thisInd = in1not2[i];
            p_obs[thisInd] += d_p1;
            
        }
        
        for(int i = 0;  i < k2; i++){
            thisInd = in2not1[i];
            p_obs[thisInd] += d_p2;
            
        }
        p_mass[p1_ind] = new_p1;
        p_mass[p2_ind] = new_p2;
        return;
    }
    
    if(k2 == 0){
        double new_p1 = p1 + p2;
        double new_p2 = 0;
        
        double d_p1 = p2;
        double d_p2 = -p2;
        for(int i = 0; i < k1; i++){
            thisInd = in1not2[i];
            p_obs[thisInd] += d_p1;
            
        }
        
        for(int i = 0;  i < k2; i++){
            thisInd = in2not1[i];
            p_obs[thisInd] += d_p2;
            
        }
        p_mass[p1_ind] = new_p1;
        p_mass[p2_ind] = new_p2;
        return;
        return;
    }   */
    
    if(k1 > 0 && k2 > 0){
    
        double S1 = 0;
        double S2 = 0;
    
        for(int i = 0; i < k1; i++){
            thisInd = in1not2[i];
            b1 = min(b1, p_obs[thisInd]);
            S1 += 1/p_obs[thisInd];
        }
        b1 -= p1;
        if(b1 == 2.0){b1 = 0.0;}
        S1 *= (b1 + p1);
    
        for(int i = 0; i < k2; i++){
            thisInd = in2not1[i];
            b2 = min(b2, p_obs[thisInd]);
            S2 += 1 / p_obs[thisInd];
        
        }
        b2 -= p2;
        if(b2 == 2.0){b2 = 0.0;}
        S2 *= (b2 + p2);
        if(S1 == 0 && S2 == 0) {return;}
    
    new_p1 = max(0.0, min(b0, (b0 + b1 + b2) * S1/(S1 + S2) - b1 ) );
    new_p2 = b0 - new_p1;
   
        
        
        if(ISNAN(new_p1) || ISNAN(new_p2)){
            Rprintf("Warning: ISNAN new_p1, new_p2! p1 = %f, p2 = %f, S1 = %f, S2 = %f, b0 = %f, ",
                    p1, p2, S1, S2, b0);
            Rprintf("b1 = %f, b2 = %f\n", b1, b2);
        }
    
    }
    else{
        if(k1 == 0){
            new_p1 = 0;
            new_p2 = p1 + p2;
        }
        else{
            new_p1 = p1 + p2;
            new_p2 = 0;
        }
    }
    double d_p1 = new_p1 - p1;
    
    if(d_p1 == 0) {return;}
    
    double d_p2 = new_p2 - p2;
    
    
    for(int i = 0; i < k1; i++){
        thisInd = in1not2[i];
        p_obs[thisInd] += d_p1;
        
    }

    for(int i = 0;  i < k2; i++){
        thisInd = in2not1[i];
        p_obs[thisInd] += d_p2;

    }
    p_mass[p1_ind] = new_p1;
    p_mass[p2_ind] = new_p2;
 
}


void SEXPMat2pmass_info(SEXP r_mat, std::vector< std::vector<int> > &cliqMat, bool obs_as_rows ){
    std::vector<int> dims = getSEXP_MatDims(r_mat);
    int N  = 0;
    int K  = 0;
    
    if(obs_as_rows == true){
        N = dims[0];
        K = dims[1];
    }
    else{
        N = dims[1];
        K = dims[0];
    }
    
    int* r_mat_ptr = INTEGER(r_mat);
    
    cliqMat.resize(K);
    int k1, k2, matInd;
    
    if(obs_as_rows == true){
        for(int i = 0; i < K; i++){
            k1 = 0;
            matInd = i * N;
            for(int j = 0; j < N; j++){
                if( r_mat_ptr[matInd] == 1) {k1++;}
                matInd++;
            }
            cliqMat[i].resize(k1);
            k2 = 0;
            matInd = i * N;
            for(int j = 0; j < N; j++){
                if( r_mat_ptr[matInd] == 1) {
                    cliqMat[i][k2] = j;
                    k2++;
                }
                matInd++;
            }
        }
    }
    else{
        std::vector<int> v_sizes(K, 0);
        
        for(int i = 0; i < N; i++){
            matInd = i * K;
            for(int j = 0; j < K; j++){
                v_sizes[j] += r_mat_ptr[matInd];
                matInd++;
            }
        }
        
        for(int j = 0; j < K; j++) { cliqMat[j].reserve(v_sizes[j]); }
        for(int i = 0; i < N; i++){
            matInd = i*K;
            for(int j = 0; j < K; j++){
                if(r_mat_ptr[matInd] == 1){
                    cliqMat[j].push_back(i);
                }
                matInd++;
    
            }
        }
    }
}



void bvcen::drop_pmass(int p_ind){
    if(p_mass[p_ind] > 0){
        Rprintf("warning: attempting to drop active point with positive mass\n");
        return;
    }
/*    int n = pmass_in_ob.size();
    for(int i = 0; i < n; i++){ drop_index(p_ind, pmass_in_ob[i]); }    */
    drop_index(p_ind, pos_pmass);
    add_index(p_ind, zer_pmass);
    int new_k = pos_pmass.size();
    dp_act.resize(new_k);
}

void bvcen::add_pmass(int p_ind){
/*    int n = pmass_in_ob.size();
    for(int i = 0; i < n; i++){ add_index(p_ind, pmass_in_ob[i]); } */
    if(p_mass[p_ind] > 0){
        Rprintf("warning: attempting to add active point, but it already has positive mass\n");
        return;
    }
    add_index(p_ind, pos_pmass);
    drop_index(p_ind, zer_pmass);
    int new_k = pos_pmass.size();
    dp_act.resize(new_k);
}

void bvcen::drop_zeros(){
    int num_pos = pos_pmass.size();
    drop_inds.resize(0);
    int ind;
    for(int i = 0; i < num_pos; i++){
        ind = pos_pmass[i];
        if(p_mass[ind] == 0) {drop_inds.push_back(ind);}
    }
    int num_drops = drop_inds.size();
    for(int i = 0; i < num_drops; i++){
        drop_pmass(drop_inds[i]);
    }
}

void bvcen::add_points(){
    int num_zeros = zer_pmass.size();
    drop_inds.resize(0);
    double cut = 1.0;
    
    calc_full_dp();
    int ind;
    for(int i = 0; i < num_zeros; i++){
        ind = zer_pmass[i];
        if(dp_full[ind] > cut){drop_inds.push_back(ind);}
    }
    int num_adds = drop_inds.size();
    for(int i = 0; i < num_adds; i++){
        add_pmass(drop_inds[i]);
    }
    vem_act();
    fullError = actError;
}

double bvcen::llk(){
    double ans = 0;
    int n = p_obs.size();
    for(int i = 0; i < n; i++){ ans += log(p_obs[i]);}
    return(ans);
}

void bvcen::update_pobs(){
    int n = p_obs.size();
    int num_pos = pos_pmass.size();
    int this_nobs, this_p_ind;
    double this_p;
    std::vector<int>* indVec;
    for(int i = 0; i < n; i++) { p_obs[i] = 0;}
    for(int i = 0; i < num_pos; i++){
        this_p_ind = pos_pmass[i];
        this_p = p_mass[this_p_ind];
        
        indVec = &(pmass_in_ob[this_p_ind]);
        this_nobs = (*indVec).size();
        for(int j = 0; j < this_nobs; j++){
            p_obs[ (*indVec)[j] ] += this_p;
        }
    }
}


void bvcen::full_em(){
    calc_full_dp();
    int n_mi = dp_full.size();
//    int n_obs = p_obs.size();
    fullError = 0;
    
    for(int i = 0; i < n_mi; i++){
        p_mass[i] *= dp_full[i];
        fullError = max(fullError, dp_full[i]);
        }
    fullError -= 1;
    update_pobs();
}

void bvcen::act_em(){
    calc_act_dp();
    int n_mi = dp_act.size();
  //  int n_obs = p_obs.size();
    actError = 0;
    int ind;
    for(int i = 0; i < n_mi; i++){
        ind = pos_pmass[i];
        p_mass[ind] *= dp_act[i];
        actError = max(actError, dp_act[i]);
    }
    actError -= 1;
    update_pobs();
}


void bvcen::calc_full_dp(){
    int n_mi = dp_full.size();
    int n_obs = p_obs.size();
    double inv_nobs = 1.0/n_obs;
    
    p_obs_inv.resize(n_obs);
    for(int i = 0; i < n_obs; i++){ p_obs_inv[i] = 1/p_obs[i]; }
    std::vector<int>* indVec;
    int this_nobs;
    int ind;
    double this_d;
    for(int i = 0; i < n_mi; i++){
        indVec = &(pmass_in_ob[i]);
        this_nobs = (*indVec).size();
        this_d = 0;
        for(int j = 0; j < this_nobs; j++){
            ind = (*indVec)[j];
            this_d += p_obs_inv[ind];
        }
        dp_full[i] = this_d;
        dp_full[i] *= inv_nobs;
    }
}


void bvcen::calc_act_dp(){
    int n_mi = pos_pmass.size();
    int n_obs = p_obs.size();
    double inv_nobs = 1.0/n_obs;
    p_obs_inv.resize(n_obs);
    for(int i = 0; i < n_obs; i++){ p_obs_inv[i] = 1/p_obs[i]; }
    std::vector<int>* indVec;
    int this_nobs;
    int this_mi;
    int ind;
    double thisError = -1;
    for(int i = 0; i < n_mi; i++){
        this_mi = pos_pmass[i];
        indVec = &(pmass_in_ob[this_mi]);
        this_nobs = (*indVec).size();
        dp_act[i] = 0;
        for(int j = 0; j < this_nobs; j++){
            ind = (*indVec)[j];
            dp_act[i] += p_obs_inv[ind];
            //SHOULD be dp_act[i], NOT dp_act[ind]
        }
        dp_act[i] *= inv_nobs;
        thisError = max(thisError, dp_act[i]);
    }
    actError = thisError - 1.0;
}

void bvcen::vem_act(){
    calc_act_dp();
    double cut = 1.0;
    int* min = new int;
    int* max = new int;
    getRelValIndices(cut, dp_act, pos_pmass, posInds, negInds, max, min);

    if( (*min) >= 0 && (*max) >= 0){ squeeze(*min, *max); }
    
    delete min;
    delete max;
    
    int i2 = -1;
    int ind1, ind2;
    int k1 = posInds.size();
    int k2 = negInds.size();
    
    if(k1 == 0 || k2 == 0){return;}
    
    for(int i = 0; i < k1; i++){
        i2++;
        if(i2 == k2) {i2 = 0;}
        ind1 = posInds[i];
        ind2 = negInds[i2];
        squeeze(ind1, ind2);
    }
    i2 = -1;
    for(int i = 0; i < k2; i++){
        i2++;
        if(i2 == k1){i2 = 0;}
        ind1 = negInds[i];
        ind2 = posInds[i2];
        squeeze(ind1, ind2);
    }
    
    int num_act = pos_pmass.size();
    for(int i = 1; i < num_act; i++){
        ind1 = pos_pmass[i-1];
        ind2 = pos_pmass[i];
        squeeze(ind1, ind2);
    }
}


bvcen::bvcen(SEXP cliqMat, SEXP R_obs_as_rows){
    std::vector<int> dims = getSEXP_MatDims(cliqMat);
    bool obs_as_rows = LOGICAL(R_obs_as_rows)[0] == TRUE;
    int n_obs = 0;
    int n_mi = 0;
    if(obs_as_rows == true){
        n_obs = dims[0];
        n_mi = dims[1];
    }
    else{
        n_obs = dims[1];
        n_mi = dims[0];
    }
    p_obs.resize(n_obs);
    SEXPMat2pmass_info(cliqMat, pmass_in_ob, obs_as_rows);
    pos_pmass.resize(n_mi);
    for(int i = 0; i < n_mi; i++){pos_pmass[i] = i;}
    p_mass.resize(n_mi, (1.0/n_mi) );
    dp_full.resize(n_mi);
    dp_act.resize(n_mi);
    update_pobs();
    
    zer_pmass.reserve(n_mi);
    drop_inds.reserve(n_mi);
    
    negInds.reserve(n_mi);
    posInds.reserve(n_mi);
    
    fullError = R_PosInf;
    actError = R_PosInf;
}

double bvcen::get_ptot(){
    double ptot = 0;
    int k = p_mass.size();
    for(int i = 0; i < k; i++){ptot += p_mass[i];}
    return(ptot);
}

/*      R TESTING FUNCTION          */

SEXP R_testDiffStep(SEXP in1, SEXP in2){
    std::vector<int> c_in1, c_in2, c_in1not2, c_in2not1;
    SEXPIndex2intIndex(in1, c_in1);
    SEXPIndex2intIndex(in2, c_in2);
    
    findIndexDiffs(c_in1, c_in2, c_in1not2, c_in2not1);
    
    
    return(R_NilValue);
}


    


SEXP optCliq(SEXP cliqMat, SEXP R_tol, SEXP R_innerLoops, SEXP R_outerLoops, SEXP R_obs_as_rows){

    double tol = REAL(R_tol)[0];
    bvcen bvObj(cliqMat, R_obs_as_rows);
    int innerLoops = INTEGER(R_innerLoops)[0];
    int outLoops = INTEGER(R_outerLoops)[0];
    int in_it = 0;
    int out_it = 0;
    int tot_its = 0;
    
    while(out_it < outLoops && bvObj.fullError > tol){
        out_it++;
        in_it = 0;
        while(in_it <innerLoops && bvObj.actError > tol){
            tot_its++;
            in_it++;
            bvObj.act_em();
            bvObj.vem_act();
            bvObj.drop_zeros();
        }
        bvObj.add_points();
    }
    
    int n_mi = bvObj.p_mass.size();
    
    SEXP pvec = PROTECT(Rf_allocVector(REALSXP, n_mi));
    for(int i = 0; i < n_mi; i++){
        REAL(pvec)[i] = bvObj.p_mass[i];
    }
    
    SEXP llh = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(llh)[0] = bvObj.llk();
    
    SEXP R_tot_its = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(R_tot_its)[0] = tot_its;
    
    SEXP out_its = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(out_its)[0] = out_it;
    
    SEXP error = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(error)[0] = bvObj.fullError;
    
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, pvec);
    SET_VECTOR_ELT(ans, 1, llh);
    SET_VECTOR_ELT(ans, 2, R_tot_its);
    SET_VECTOR_ELT(ans, 3, out_its);
    SET_VECTOR_ELT(ans, 4, error);
    UNPROTECT(6);
    
    return(ans);
}
