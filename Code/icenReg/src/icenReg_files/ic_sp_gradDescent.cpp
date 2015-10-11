//
//  ic_sp_gradDescent.cpp
//  
//
//  Created by Cliff Anderson Bergman on 10/9/15.
//
//

/*   Gradient Ascent Step */

void icm_Abst::baseCH_2_baseS(){
    int k = baseCH.size();
    baseS.resize(k);
    baseS[0] = 1.0;
    baseS[k-1] = 0.0;
    for(int i = 1; i < (k-1); i++){
        baseS[i] = exp(-exp(baseCH[i]));
    }
}

void icm_Abst::baseS_2_baseP(){
    int k = baseS.size() - 1;
    baseP.resize(k);
    for(int i = 0; i < k; i++){
        baseP[i] = baseS[i] - baseS[i+1];
    }
}

void icm_Abst::baseP_2_baseS(){
    int k = baseP.size();
    baseS.resize(k + 1);
    baseS[0] = 1.0;
    for(int i = 1; i < k; i++){
        baseS[i] = baseS[i-1] - baseP[i-1];
    }
    baseS[k] = 0.0;
}

void icm_Abst::baseS_2_baseCH(){
    int k = baseS.size();
    baseCH.resize(k);
    baseCH[0] = R_NegInf;
    baseCH[k-1] = R_PosInf;
    for(int i = 1; i < (k-1); i++){
        baseCH[i] = log(-log(baseS[i]));
    }
}

double icm_Abst::llk_from_p(){
    baseP_2_baseS();
    baseS_2_baseCH();
    double ans = sum_llk();
    return(ans);
}


double icm_Abst::dervConS_fromBaseS(double s, double eta){
    bool adjust_h = false;
    double org_h = h;
    if(h > s*2 || h > (1-s)*2 ){
        h = s/2.0;
        if(s > 0.5){
            h = (1.0-s)/2.0;
        }
        adjust_h = true;
    }
    double f_h, f_l,  ans;
    s += h;
    f_h = baseS2CondS(s, eta);
    s -= 2*h;
    f_l = baseS2CondS(s, eta);
    ans = (f_h - f_l) / (2 * h);
    if(adjust_h){
        h = org_h;
    }
    if(ISNAN(ans) || ans == R_NegInf || ans == R_PosInf){
        ans = 0;
    }
    return(ans);
}

void icm_Abst::calc_cond_S_derv(){
    double this_eta, this_p_obs;
    double s_l, s_r;
    int n = expEtas.size();
    int k = baseCH.size();
    int l_ind, r_ind;
    d_cond_S_left.resize(k);
    d_cond_S_right.resize(k);
    for(int i = 0; i < k; i++){
        d_cond_S_left[i] = 0;
        d_cond_S_right[i] = 0;
    }
    for(int i = 0; i < n; i++){
        l_ind = obs_inf[i].l;
        r_ind = obs_inf[i].r;
        this_eta = etas[i];
        this_p_obs = obs_inf[i].pob;
        if(l_ind > 0){
            s_l = exp( - exp( baseCH[l_ind]) );
            d_cond_S_left[l_ind - 1] -= dervConS_fromBaseS( s_l , this_eta) / this_p_obs;
        }
        if(r_ind < (k-1) ){
            s_r = exp( - exp( baseCH[r_ind + 1]) );
            d_cond_S_right[r_ind] += dervConS_fromBaseS( s_r , this_eta) / this_p_obs;
        }
    }
}


/*  void icm_Abst::makeConS_Sum(){
    int k = baseCH.size();
    base_p_derv.resize(k - 1);
    base_p_derv[0] = d_cond_S_left[0] + d_cond_S_right[0];
    for(int i = 1; i < (k-1); i++){
        base_p_derv[i] = base_p_derv[i - 1] + d_cond_S_left[i] + d_cond_S_right[i];
    }
}   */

void icm_Abst::calc_base_p_derv(){
    calc_cond_S_derv();
    int k = baseCH.size() - 1;
    base_p_derv.resize(k);
    base_p_derv[ k-1 ] = 0;
    int k_l, k_r, l_ind, r_ind;
    k_l = node_inf[ k-1 ].l.size();
    k_r = node_inf[ k-1 ].r.size();
    for(int j = 0; j < k_l; j++){
        l_ind = node_inf[ k-1 ].l[j];
        base_p_derv[ k-1 ] += d_cond_S_left[l_ind];
    }
    for(int j = 0; j < k_r; j++){
        r_ind = node_inf[ k-1 ].r[j];
        base_p_derv[ k-1 ] += d_cond_S_right[r_ind + 1];
    }
    
    for(int i = k - 2; i >= 0; i--){
        base_p_derv[i] = base_p_derv[i+1];
        k_l = node_inf[i].l.size();
        k_r = node_inf[i].r.size();
        
        for(int j = 0; j < k_l; j++){
            l_ind = node_inf[i].l[j];
            base_p_derv[i] += d_cond_S_left[l_ind];
        }
        for(int j = 0; j < k_r; j++){
            r_ind = node_inf[i].r[j];
            base_p_derv[i] += d_cond_S_right[r_ind + 1];
        }

        
    }
    
    Rprintf("3 base derivative = %f, %f, %f\n", base_p_derv[0], base_p_derv[1], base_p_derv[base_p_derv.size()-1]);
}

double icm_Abst::getMaxScaleSize(vector<double> &p, vector<double> &prop_p){
    double max_scale = 2.0;
    int k = p.size();
    int k2 = prop_p.size();
    if(k != k2){
        Rprintf("warning: k != k2 in getMaxScaleSize k = %d, k2 = %d\n", k, k2);
        return(0.0);
    }
    double thisMax = 1.0;
    for(int i = 0; i < k; i++){
        if(prop_p[i] != 0 && p[i] > 0){
            thisMax = max(-p[i]/prop_p[i], (1.0 - p[i]) / prop_p[i]);
            if(ISNAN(thisMax)){
                Rprintf("\nWarning: ISNAN(thisMax). p[i] = %f, prop_p[i] = %f\n",
                        p[i], prop_p[i]);
            }
        }
        max_scale = min(max_scale, thisMax);
    }
    
    return(max_scale);
}


void icm_Abst::gradientDescent_step(){
    baseCH_2_baseS();
    baseS_2_baseP();
    
    calc_base_p_derv();
    int k = base_p_derv.size();
    prop_p.resize(k);
    double prop_mean = 0;
    int act_sum = 0;
    
    vector<bool> isActive(k);
    for(int i = 0; i < k; i++){
        if(base_p_derv[i] >= 0 || baseP[i] > 0){
            isActive[i] = true;
            act_sum++;
        }
        else { isActive[i] = false; }
    }
    
    for(int i = 0; i < k; i++){
        if(isActive[i]){ prop_mean += base_p_derv[i]; }
    }
    
    prop_mean = prop_mean / act_sum;
    
    for(int i = 0; i < k; i++){
        if(isActive[i]){ prop_p[i] = base_p_derv[i] - prop_mean;}
        else {prop_p[i] = 0.0;}
    }
    
    
    
    double scale_max = getMaxScaleSize(baseP, prop_p);
    double delta_val = scale_max/2.0;
    
    
    if(delta_val == 0){
        Rprintf("delta_val = 0, quitting gradientDescent_step\n");
        return;
    }
    
    
    
    mult_vec( delta_val, prop_p );
    add_vec(prop_p, baseP);
    double llk_hh = llk_from_p();
    mult_vec( -1.0/2.0, prop_p );
    add_vec(prop_p, baseP);
    double llk_h = llk_from_p();
    add_vec(prop_p, baseP);
    double llk_0 = llk_from_p();
    
    for(int i = 0; i < k; i++){
        if(isActive[i]){ prop_p[i] = base_p_derv[i] - prop_mean; }
        else{ prop_p[i] = 0.0;}
    }
    
    double d1 = ( llk_hh - llk_0 ) / ( delta_val );
    double d2 = (llk_hh + llk_0 - 2.0 * llk_h ) / (delta_val * delta_val / 4);
    delta_val = -d1/d2;
    
    double analyticDircDerv = directional_derv(base_p_derv, prop_p);
    Rprintf("Direct Derivative: Analytic: %f, Numeric: %f\n", analyticDircDerv, d1);
    
    if(delta_val <= 0){
        Rprintf("note: delta_val <= 0, equal to %f. d1 = %f, Quitting GD step.\n", delta_val, d1);
        return;
    }
    
    delta_val = min( delta_val, scale_max );
    
    mult_vec(delta_val, prop_p);
    add_vec(prop_p, baseP);
    
    double new_llk = llk_from_p();
    mult_vec(-1.0, prop_p);
    int tries = 0;
    while(tries < 50 && new_llk  < llk_0){
        tries++;
        mult_vec(0.5, prop_p);
        add_vec(prop_p, baseP);
        new_llk = llk_from_p();
    }
    if(new_llk < llk_0){
        add_vec(prop_p, baseP);
        new_llk = llk_from_p();
        if(new_llk < llk_0){
            Rprintf("warning: new_llk < llk_0 in the very end of GD step! Diff = \n",
                    new_llk - llk_0);
        }
    }
}

