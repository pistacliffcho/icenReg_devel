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
                thisMax = 1.0;
            }
        }
        max_scale = min(max_scale, thisMax);
    }
    
    return(max_scale);
}

void icm_Abst::EM_step(){
	double org_llk = sum_llk();
	
    backupCH = baseCH;
    baseCH_2_baseS();
    baseS_2_baseP();

    numeric_dobs_dp(false);
    int k = base_p_derv.size();
//    double n = etas.size();
	baseP_backup.resize(k);
	for(int i = 0; i < k; i++){
		baseP_backup[i] = baseP[i];
//		baseP[i] *= (n + base_p_derv[i]);
		baseP[i] *= (base_p_derv[i]);
		if(baseP[i] < 0){baseP[i] = 0;}
	}
	double sum_p = 0;
	for(int i = 0; i < k; i++){ sum_p += baseP[i]; }
	for(int i = 0; i < k; i++){ baseP[i] /= sum_p; }	
	
	double new_llk = llk_from_p();	
	
	if(new_llk < org_llk){
	//	Rprintf("Note: EM step failed. Difference in llk = %f, current iter = %d\n", new_llk - org_llk, iter);
		for(int i = 0; i < k; i++){ baseP[i] = baseP_backup[i];}
		new_llk = llk_from_p();
	}
}

void icm_Abst::numeric_dobs2_d2p(){
        	
    backupCH = baseCH;
    baseCH_2_baseS();
    baseS_2_baseP();
    
    double offSet = h * 0.00001;
    int k = baseP.size();
    for(int i = 0; i < k; i++){ baseP[i]+= offSet; }
	baseP_2_baseS();
    	
    numeric_dobs_dp(true);
    k = base_p_derv.size();
    base_p_derv2.resize(k);
    base_p_2ndDerv.resize(k);
    for(int i = 0; i < k; i++){ base_p_derv2[i] = base_p_derv[i]; }
    for(int i = 0; i < k; i++){ baseP[i]-= 2.0 * offSet; }
    baseP_2_baseS();
    numeric_dobs_dp(true);
	for(int i = 0; i < k; i++){ base_p_2ndDerv[i] = (base_p_derv2[i] - base_p_derv[i]) / (2.0 * offSet); }
	for(int i = 0; i < k; i++){ base_p_derv[i] = (base_p_derv2[i] + base_p_derv[i]) / 2.0; }


    for(int i = 0; i < k; i++){ baseP[i]+= offSet; }
	baseP_2_baseS();
	
}

void icm_Abst::experimental_step(){
    
	if(failedGA_counts > 500){return;}
	
	double org_llk = sum_llk();
	
    backupCH = baseCH;
    baseCH_2_baseS();
    baseS_2_baseP();
    	
    numeric_dobs2_d2p();
    int k = base_p_derv.size();
    prop_p.resize(k);
    double prop_mean = 0;
    int act_sum = 0;
    double new_llk;
    
    vector<bool> isActive(k);
    for(int i = 0; i < k; i++){
        if(baseP[i] > 0 && !ISNAN(base_p_derv[i]) && base_p_2ndDerv[i] < -0.001){
            isActive[i] = true;
            act_sum++;
        }
        else { isActive[i] = false; }
    }
    
    for(int i = 0; i < k; i++){
        if(isActive[i]){ prop_mean += -base_p_derv[i]/base_p_2ndDerv[i]; }
    }
    
    
    prop_mean = prop_mean / act_sum;
    
    for(int i = 0; i < k; i++){
        if(isActive[i]){ prop_p[i] = -base_p_derv[i]/base_p_2ndDerv[i] - prop_mean;}
        else {prop_p[i] = 0.0;}
    }
    
    
    makeUnitVector(prop_p);
    
    
    double scale_max = getMaxScaleSize(baseP, prop_p);

    
    for(int i = 0; i < k; i++){ prop_p[i] *= -1.0; }
    scale_max = min(scale_max, getMaxScaleSize(baseP, prop_p));
    for(int i = 0; i < k; i++){ prop_p[i] *= -1.0; }
    
    double delta_val = scale_max/2.0;
    
    delta_val = min(delta_val, h);
    delta_val = delta_val/10.0;
    
//    double analytic_dd = directional_derv(base_p_derv, prop_p);
    
    
    if(delta_val == 0){
        failedGA_counts++;
        baseCH = backupCH;
        new_llk = sum_llk();
        Rprintf("Exit 1\n");
        return;
    }
    
    add_vec(delta_val, prop_p, baseP);
    double llk_h = llk_from_p();
    add_vec(-2.0 * delta_val, prop_p, baseP);
    double llk_l = llk_from_p();
    add_vec(delta_val, prop_p, baseP);
    double llk_0 = llk_from_p();
    
	
    double d1 = ( llk_h - llk_l ) / ( 2 * delta_val );
    double d2 = (llk_h + llk_l - 2.0 * llk_0 ) / (delta_val * delta_val);
        
    delta_val = -d1/d2;
	    
    if(ISNAN(delta_val)){
        failedGA_counts++;
        baseCH= backupCH;
        new_llk = sum_llk();
        Rprintf("warning: delta_val is nan in GA step. llk_h = %f, llk_l = %f, llk_0 = %f, scale_max = %f\n", 
    			llk_h, llk_l, llk_0, scale_max);
        Rprintf("Exit 3\n");
        return;
    }
    
    scale_max = getMaxScaleSize(baseP, prop_p);
    delta_val = min( delta_val, scale_max );
    
    add_vec(delta_val, prop_p, baseP);

    new_llk = llk_from_p();
    mult_vec(-1.0, prop_p);
    int tries = 0;
    
    double this_delta = delta_val;
    
    while(tries < 5 && new_llk  < llk_0){
        tries++;
        this_delta = this_delta/2;
        add_vec(this_delta, prop_p, baseP);
        new_llk = llk_from_p();
    }
    if(new_llk < llk_0){
        failedGA_counts++;
        baseCH = backupCH;
        new_llk = sum_llk(); //Should NOT be llk_from_p(), since we are resetting the CH
        Rprintf("Exit 4\n");
		return;
    }
	
	if(org_llk > new_llk){
		failedGA_counts++;
		baseCH = backupCH;
		new_llk = sum_llk();
	}
		
}


void icm_Abst::gradientDescent_step(){
    
	if(failedGA_counts > 500){return;}
	
	double org_llk = sum_llk();
	
    backupCH = baseCH;
    baseCH_2_baseS();
    baseS_2_baseP();
    	
    numeric_dobs_dp(true);
    int k = base_p_derv.size();
    prop_p.resize(k);
    double prop_mean = 0;
    int act_sum = 0;
    double new_llk;
    
    vector<bool> isActive(k);
    for(int i = 0; i < k; i++){
        if(baseP[i] > 0 && !ISNAN(base_p_derv[i]) ){
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
    
    
    makeUnitVector(prop_p);
    
    double scale_max = getMaxScaleSize(baseP, prop_p);

    
    for(int i = 0; i < k; i++){ prop_p[i] *= -1.0; }
    scale_max = min(scale_max, getMaxScaleSize(baseP, prop_p));
    for(int i = 0; i < k; i++){ prop_p[i] *= -1.0; }
    
    double delta_val = scale_max/2.0;
    
    delta_val = min(delta_val, h);
    delta_val = delta_val/10.0;
    
    double analytic_dd = directional_derv(base_p_derv, prop_p);
    
    
    if(delta_val == 0){
        failedGA_counts++;
        baseCH = backupCH;
        new_llk = sum_llk();
        return;
    }
    
    add_vec(delta_val, prop_p, baseP);
    double llk_h = llk_from_p();
    add_vec(-2.0 * delta_val, prop_p, baseP);
    double llk_l = llk_from_p();
    add_vec(delta_val, prop_p, baseP);
    double llk_0 = llk_from_p();
    
	
    double d1 = ( llk_h - llk_l ) / ( 2 * delta_val );
    double d2 = (llk_h + llk_l - 2.0 * llk_0 ) / (delta_val * delta_val);
    
    if(iter % 2 ==0){ d1 = analytic_dd; }
    
    delta_val = -d1/d2;
	
    if(!(delta_val > 0)){
        failedGA_counts++;
        baseCH = backupCH;
        new_llk = sum_llk();
        return;
    }

    
    if(ISNAN(delta_val)){
        failedGA_counts++;
        baseCH= backupCH;
        new_llk = sum_llk();
        return;
    }

    
    scale_max = getMaxScaleSize(baseP, prop_p);
    delta_val = min( delta_val, scale_max );
    
    add_vec(delta_val, prop_p, baseP);

    new_llk = llk_from_p();
    mult_vec(-1.0, prop_p);
    int tries = 0;
    
    double this_delta = delta_val;
    
    while(tries < 5 && new_llk  < llk_0){
        tries++;
        this_delta = this_delta/2;
        add_vec(this_delta, prop_p, baseP);
        new_llk = llk_from_p();
    }
    if(new_llk < llk_0){
        failedGA_counts++;
        baseCH = backupCH;
        new_llk = sum_llk(); //Should NOT be llk_from_p(), since we are resetting the CH
		return;
    }
	
	if(org_llk > new_llk){
		failedGA_counts++;
		baseCH = backupCH;
		new_llk = sum_llk();
	}
		
//	Rprintf("change in llk in CGA step = %f\n", new_llk - org_llk);	
		
}




double icm_Abst::cal_log_obs(double s1, double s2, double eta){
    double l = baseS2CondS(s1, eta);
    double r = baseS2CondS(s2, eta);
    return(log(l - r) );
}



void icm_Abst::numeric_dobs_dp(bool forGA){    
    int p_k = baseS.size();
    int k = p_k - 1;
    int n = etas.size();
    dob_dp_both.resize(n);
    dob_dp_rightOnly.resize(n);
    int lind, rind;
	double h_mult = 0.0001;
   	h *= h_mult;

	if(forGA){
	    double sl, sr, llk_h,llk_l, this_eta, this_h;    
   	 
   		for(int i = 0; i < n; i++){
    	    sl = baseS[ obs_inf[i].l];
    	    sr = baseS[ obs_inf[i].r + 1];
    	    this_eta = etas[i];
    	    if(sl == 1.0 && sr == 0.0){
    	        dob_dp_rightOnly[i] = 0;
    	        dob_dp_both[i] = 0;
    	    }
    	    else if(sr == 0){
    	        dob_dp_rightOnly[i] = 0;
    	        this_h = min(sl/2.0, h);
    	        sl -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sl += this_h * 2.0;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[i] = (llk_h - llk_l) / (2 * this_h);
    	    }
    	    else if( sl == 1.0 ){
    	        this_h = min(sr / 2.0, h);
    	        sr -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sr += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[i] = (llk_h - llk_l)/(2*this_h);
    	        dob_dp_rightOnly[i] = dob_dp_both[i];
    	    }
    	    else{
    	        this_h = min(sr /2.0, h);
    	        sr -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sr += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        sr -= this_h;
    	        dob_dp_rightOnly[i] = (llk_h - llk_l)/(2*this_h);
    	        sr -= this_h;
    	        sl -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
            
    	        sr += 2.0 * this_h;
    	        sl += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[i] = (llk_h - llk_l)/(2*this_h);
    	        
    	    }
    	}
    }
    else{
    	for(int i = 0; i < p_k; i++){
    		dob_dp_both[i] = 0;
    		dob_dp_rightOnly[i] = 0;
    	}
    	double thisProb;
    	double num_n = n;
    	for(int i = 0; i < n; i++){
    		lind = obs_inf[i].l;
    		rind = obs_inf[i].r + 1;
    		thisProb = baseS[lind] - baseS[rind];
    		dob_dp_rightOnly[i] = 1.0/(num_n*thisProb);
    	}
    }
    
    base_p_derv.resize(k);
	
    int k_l, k_r;
    node_info* nd;
    for(int j = k-1; j >=0; j--){
        nd = &node_inf[j+1];
        k_r = nd->r.size();
        k_l = nd->l.size();
        if(j != k-1){
            base_p_derv[j] = base_p_derv[j+1];
        }
        else{
            base_p_derv[j] = 0;
        }
        for(int i = 0; i < k_r; i++){
            rind = nd->r[i];
            base_p_derv[j] += dob_dp_rightOnly[rind] * w[rind];
        }
        for(int i = 0; i < k_l; i++){
            lind = nd->l[i];
            base_p_derv[j] -= dob_dp_rightOnly[lind] * w[lind];
            base_p_derv[j] += dob_dp_both[lind] * w[lind];
        }
    }

    h = h/h_mult;

}

