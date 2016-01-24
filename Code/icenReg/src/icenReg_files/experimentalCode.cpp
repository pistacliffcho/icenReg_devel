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

    numeric_dobs_dp(true);
    
    int min_ind, max_ind;
    min_ind = 0;
    max_ind = 0;
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
    
    exchange_p_opt(max_ind, min_ind);
    
}


void icm_Abst::vem_sweep(){
    
    baseCH_2_baseS();
    baseS_2_baseP();
    
    numeric_dobs_dp(true);
    int k = baseP.size();
    int ind1, ind2;
    bool foundPos = false;
    double meanDervs = 0;
    double sumAct = 0;
    for(int i = 0; i < (k-1); i++){
        if(baseP[i] > 0){
            sumAct++;
            meanDervs += base_p_derv[i];
        }
    }
    
    meanDervs = meanDervs/sumAct;
    for(int i = 0; i < k; i++){
        if(!foundPos){
            if(base_p_derv[i] > 0 && baseP[i] > 0){
                foundPos = true;
                ind1 = i;
            }
        }
        else{
            if(base_p_derv[i] < 0 && baseP[i] > 0){
                foundPos = false;
                ind2 = i;
                exchange_p_opt(ind1 , ind2);
            }
        }
    }
}

void icm_Abst::vem_sweep2(){
    
    baseCH_2_baseS();
    baseS_2_baseP();
    int ind1, ind2;
    bool found1 = false;
    int k = baseP.size();
    for(int i = 0; i < k; i++){
        if(!found1){
            if(baseP[i] >0){
                found1 = true;
                ind1 = i;
            }
        }
        else{
            if(baseP[i] > 0){
                ind2 = i;
                exchange_p_opt(ind1 , ind2);
                if(baseP[ind2] <= 0){
                    found1 = false;
                }
                else{
                    ind1 = ind2;
                }
            }
        }
    }
}




void icm_Abst::exchange_p_opt(int i1, int i2){
    double minLoss = min(baseP[i1] , baseP[i2]);
    double this_h = min(h, minLoss / 10.0);
    
    if(!(this_h > 0)){
        return;
    }
    
    double llk_h = exchangeAndUpdate(this_h, i1, i2);
    double llk_l = exchangeAndUpdate(-2.0 * this_h, i1, i2);
    double llk_0 = exchangeAndUpdate(this_h, i1, i2);
    
    
    double d1 = (llk_h - llk_l) / (2 * this_h);
    double d2 = (llk_h + llk_l - 2.0*llk_0) / (this_h * this_h);
    
    double delta = -d1/d2;
    
    if(delta < -baseP[i1]){
        delta = -baseP[i1];
    }
    if(delta > baseP[i2]){
        delta = baseP[i2];
    }
    
    if(!(d2 < 0) || ISNAN(delta) || delta == R_PosInf || delta == R_NegInf ){
        return;
    }
    double llk_new = exchangeAndUpdate(delta, i1, i2);
    if(llk_new < llk_0){
        llk_new = exchangeAndUpdate(-0.5 * delta, i1, i2);
        if(llk_new < llk_0){
            llk_new = exchangeAndUpdate(-0.5 * delta, i1, i2);
        }
    }
}


void getUniqInts(int p_i1, int p_i2, vector<int> &uniqInts,
                 vector<node_info> &vec_vec, vector<bool> &usedVec){
    uniqInts.clear();
    int i1 = p_i1 + 1;
    int i2 = p_i2 + 1;
    int tot_pos_lng = 0;
    int lv = vec_vec.size();
    if(i1 >= lv){Rprintf("i1 too big in getUniqInts\n");return;}
    if(i2 >= lv){Rprintf("i2 too big in getUniqInts\n");return;}
    
    for(int i = i1; i < i2; i++){
        tot_pos_lng += vec_vec[i].l.size();
        tot_pos_lng += vec_vec[i].r.size();
    }
    
    uniqInts.reserve(tot_pos_lng);
    vector<int>* theseIndices;
    int thisSize;
    int thisIndex;
    for(int i = i1; i < i2; i++){
        theseIndices = &vec_vec[i].l;
        thisSize = theseIndices->size();
        for(int j = 0; j < thisSize; j++){
            thisIndex = (*theseIndices)[j];
            if(usedVec[thisIndex] == false){
                usedVec[thisIndex] = true;
                uniqInts.push_back(thisIndex);
            }
        }
        theseIndices = &vec_vec[i].r;
        thisSize = theseIndices->size();
        for(int j = 0; j < thisSize; j++){
            thisIndex = (*theseIndices)[j];
            if(usedVec[thisIndex] == false){
                usedVec[thisIndex] = true;
                uniqInts.push_back(thisIndex);
            }
        }
    }
    thisSize = usedVec.size();
    int thatSize = uniqInts.size();
    for(int i = 0; i < thatSize; i++){
        thisIndex = uniqInts[i];
        if(thisIndex >= thisSize){
            Rprintf("warning: thisIndex >= thisSize. thisIndex = %d, thisSize = %d\n", thisIndex, thisSize);
            return;
        }
        usedVec[thisIndex] = false;
    }
}


double icm_Abst::exchangeAndUpdate(double delta, int i1, int i2){
    exchange(delta, i1, i2, baseP);
    getUniqInts(i1, i2, exchangeIndices, node_inf, usedVec);
    int thisSize = baseS.size();
    if(thisSize <= i2){
        Rprintf("warning: thisSize <= i2\n");
        return(0.0);
    }
    thisSize  = baseCH.size();
    if(thisSize <= i2){
        Rprintf("warning: thisSize <= i2-pt2\n");
        return(0.0);
    }
    for(int i = i1; i < i2; i++){
        baseS[i+1] -= delta;
        baseCH[i+1] = log(-log(baseS[i+1]));
    }
    
    double ans = 0;
    int thisInd;
    int updateSize = exchangeIndices.size();
    for(int i = 0; i < updateSize; i++){
        thisInd = exchangeIndices[i];
        update_p_ob(thisInd);
        ans += log(obs_inf[thisInd].pob) * w[thisInd];
    }
    return(ans);
}


