//
//  experimentalCode.cpp
//  
//
//  Created by Clifford Anderson-Bergman on 10/18/15.
//
//

void add_2_last(double delta, vector<double> &p){
    int k = p.size();
    double sum_others = 1.0 - p[k-1];
    double mult_fact = (sum_others - delta) / sum_others;
    
    for(int i = 0; i < (k-1); i++){
        p[i] *= mult_fact;
    }
    p[k-1] += delta;
}

void icm_Abst::last_p_update(){
    
    baseCH_2_baseS();
    baseS_2_baseP();
    
    
    double this_h = min(h, baseP[baseP.size() - 1]/10.0);
    
    add_2_last(this_h, baseP);
    double llk_h = llk_from_p();
    add_2_last(-2.0 * this_h, baseP);
    double llk_l = llk_from_p();
    add_2_last(this_h, baseP);
    double llk_0 = llk_from_p();

    double d1 = (llk_h - llk_l) / (2 * this_h);
    double d2 = (llk_h + llk_l - 2.0*llk_0) / (this_h * this_h);
    
    
    double delta = -d1/d2;
    if(!(d2 < 0) || ISNAN(delta) || delta == R_PosInf || delta == R_NegInf ){
        return;
    }
    add_2_last(delta, baseP);
    double llk_new = llk_from_p();
    if(llk_new < llk_0){
        add_2_last( -1.0 * delta, baseP );
        llk_new = llk_from_p();
    }
}

void exchange(double delta, int i, int j, vector<double> &p){
    p[i] += delta;
    p[j] -= delta;
}

void icm_Abst::vem(){

    baseCH_2_baseS();
    baseS_2_baseP();

    numeric_dobs_dp();
    int min_ind, max_ind;
    double minVal = R_PosInf;
    double maxVal = R_NegInf;
    
    int k = baseP.size();
    for(int i = 0; i < k; i++){
        if(minVal > base_p_derv[i] && baseP[i] > 0){
            minVal = base_p_derv[i];
            min_ind = i;
        }
        if(maxVal < base_p_derv[i]&& baseP[i] > 0){
            maxVal = base_p_derv[i];
            max_ind = i;
        }
    }
    
    double minLoss = baseP[min_ind];
    double this_h = min(h, minLoss / 10.0);
    this_h = min(this_h, baseP[max_ind]);
    exchange(this_h, min_ind, max_ind, baseP);
    double llk_h = llk_from_p();
    exchange(-2.0 * this_h, min_ind, max_ind, baseP);
    double llk_l = llk_from_p();
    exchange(this_h, min_ind, max_ind, baseP);
    double llk_0 = llk_from_p();

    
    double d1 = (llk_h - llk_l) / (2 * this_h);
    double d2 = (llk_h + llk_l - 2.0*llk_0) / (this_h * this_h);
    
    double delta = -d1/d2;
    if(!(d2 < 0) || ISNAN(delta) || delta == R_PosInf || delta == R_NegInf ){
        return;
    }
    exchange(delta, min_ind, max_ind, baseP);
    double llk_new = llk_from_p();
    if(llk_new < llk_0){
        exchange(-1.0 * delta, min_ind, max_ind, baseP);
        llk_new = llk_from_p();
    }
 //   Rprintf("delta llk in vem step = %f\n", llk_new - llk_0);
    
}






