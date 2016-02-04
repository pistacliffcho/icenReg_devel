#include<iostream>

/* Basic Function Call */

SEXP EMICM(SEXP Rlind, SEXP Rrind, SEXP iters){
	double old_llk = R_NegInf;
	int maxIters = 1000;
	int ems = INTEGER(iters)[0];
	double tol = pow(10.0, -10.0);
	int iter = 0;
	
	emicm emicmObj(Rlind, Rrind);
	emicmObj.llk(true);
	while( (emicmObj.current_llk - old_llk) > tol && iter < maxIters){
		iter++;
		old_llk = emicmObj.current_llk;
		emicmObj.em_step(ems);
		emicmObj.icm_step();
	}
	cout << "iters = "<< iter << "\n";
	return(R_NilValue);
}


/* EM step */

void emicm::em_step(int iters){
	p2s();
	update_p_ob(true);
	for(int i = 0; i < iters; i++){
		calc_m_for_em();
		int k = baseP.size();
		double tot = 0;
		for(int j = 0; j < k; j++){ 
			baseP[j] *= em_m[j];
			if(baseP[j] < 0) { cout << "warning: baseP[j] < 0. j = " << j << "\n";}
			tot += baseP[j];
		}
		for(int j = 0; j < k; j++){ baseP[j] /= tot; }
		p2s();
		update_p_ob(true); 
	}
}


/* ICM step */

void emicm::icm_step(){
	p2s();
	s2ch();
	calc_icm_ders();
	int k = baseCH.size() - 2; // first and last base CH are fixed at -inf and inf

	prop_delta.resize(k);
	innerCH.resize(k);
	for(int i = 0; i < k; i++){ innerCH[i] = baseCH[i+1]; }
	pavaForOptim(ch_d1, ch_d2, innerCH, prop_delta);
	
	double org_llk = llk(false);
	
	addIcmProp( baseCH, prop_delta);
	double new_llk = llk(false);

	int tries = 0;
	prop_delta *= -1.0;
	
	while(tries < 3 && new_llk < org_llk){
		tries++;
		prop_delta *= 0.5;
		addIcmProp( baseCH, prop_delta);
		new_llk = llk(false);
		
	}
	
	if(new_llk < org_llk){
		addIcmProp( baseCH, prop_delta);
		new_llk = llk(false);
	}	
	ch2p();
}




/* utilities for EM-step */


void emicm::calc_m_for_em(){
	double this_m = 0.0;
	node_info* thisNode;
	int k = baseP.size();
	double n = obs_inf.size();
	em_m.resize(k);


	thisNode = &node_inf[0];
		for(unsigned int j = 0; j < thisNode->l.size(); j++){
			this_m += 1.0 / pobs[ thisNode->l[j] ];
		} 
	em_m[0] = this_m / n;
	for(int i = 1; i < k; i++){
		thisNode = &node_inf[i];
		for(unsigned int j = 0; j < thisNode->l.size(); j++){
			this_m += 1.0 / pobs[ thisNode->l[j] ];
		} 
		thisNode = &node_inf[i-1];
		for(unsigned int j = 0; j < thisNode->r.size(); j++){
			this_m -= 1.0 / pobs[ thisNode->r[j] ];
		}
		em_m[i] = this_m / n; 
	}
}


/* utilities for ICM-step */

void icmDervParts(double lch, double* prt1, double* prt2, double* prt3){
	double ch = exp(lch);
	(*prt1) = exp(lch -ch); 
	(*prt2) = (*prt1) * (1 - ch);
	(*prt3) = (*prt1) * (*prt1);
}

double icmFirstDerv(double prt1, double pob, bool isLeft){
	if(isLeft) {return(-prt1 / pob);}
	return(prt1/pob); 
}
double icmSecondDerv(double prt2, double prt3, double pob, bool isLeft){ 
	if(isLeft) {return( -prt2 / pob - prt3 / (pob * pob) );}
	return( prt2 / pob - prt3 / (pob * pob)) ;
}

void emicm::calc_icm_ders(){
	int k = baseCH.size() - 2;
	ch_d1.resize(k);
	ch_d2.resize(k);
	
	ch_d1_con.resize(k);
	ch_d2_con1.resize(k);
	ch_d2_con2.resize(k);
	
	double thisPob;
	for(int i = 0; i < k; i++){
		icmDervParts(baseCH[i + 1], &ch_d1_con[i], 
		             &ch_d2_con1[i], &ch_d2_con2[i]);
		ch_d1[i] = 0;
		ch_d2[i] = 0;
		}
	int n = pobs.size();
	int lind, rind;
	for(int i = 0; i < n; i++){
		thisPob = pobs[i];
		lind = obs_inf[i].l;
		rind = obs_inf[i].r + 1;
		if(lind > 0){
			ch_d1[lind-1] += icmFirstDerv(ch_d1_con[lind-1], thisPob, true);
			ch_d2[lind-1] += icmSecondDerv(ch_d2_con1[lind-1], 
										   ch_d2_con2[lind-1], thisPob, true);
		}
		if(rind < (k+1) ){
			ch_d1[rind-1] += icmFirstDerv(ch_d1_con[rind-1], thisPob, false);
			ch_d2[rind-1] += icmSecondDerv(ch_d2_con1[rind-1], 
										   ch_d2_con2[rind-1], thisPob, false);
		}
	}
}

void addIcmProp(Eigen::VectorXd &bch, Eigen::VectorXd &prop){
	int k1 = bch.size();
	int k2 = prop.size();
	if(k1 != (k2 + 2) ){
		Rprintf("error: bch.size() != k2 prop.size() + 2\n");
		return;
	}
	for(int i = 1; i <= k2; i++){
		bch[i] += prop[i-1];
	}
}


void emicm::printICMdervs(){
	double llk0 = current_llk;
	double llk_l, llk_h;
	
	cout << "numeric derivatives\n";
	
	int k = prop_delta.size();
	double h = 0.0001;
	for(int i = 1; i <= k; i++){
		baseCH[i] += h;
		llk_h = llk(false);
		baseCH[i] -= 2*h;
		llk_l = llk(false);
		baseCH[i] += h;
		cout << "d1 = " << (llk_h - llk_l) / (2 * h) 
			 << "  d2 = " << (llk_h + llk_l - 2 * llk0) / (h*h) << "\n";
	}
	cout << "\n";
}


/* General utilities */

void emicm::update_p_ob(int i, bool useS){
	if(useS){ 
		pobs[i] = baseS[ obs_inf[i].l ] - baseS[ obs_inf[i].r + 1]; 
	}
	else{
		double chl = baseCH[ obs_inf[i].l ];
		double chr = baseCH[ obs_inf[i].r + 1]; 
		
		pobs[i] = exp(-exp(chl)) - exp(-exp(chr));
	}
}

void emicm::update_p_ob(bool useS){
	int n = pobs.size();
	for(int i = 0; i < n; i++){ update_p_ob(i, useS);}
}

double emicm::llk(bool useS){
	current_llk = 0;
	int n = pobs.size();
	for(int i = 0; i < n; i++){ 
		update_p_ob(i, useS); 
		current_llk += log(pobs[i]);
		}
	if(ISNAN(current_llk)){ current_llk = R_NegInf; }
	return(current_llk);
}

emicm::emicm(SEXP Rlind, SEXP Rrind){
    int n = LENGTH(Rlind);
    if(n != LENGTH(Rrind)){Rprintf("length of Rlind and Rrind not equal\n"); return;}
    
    pobs.resize(n);
    
    int* clind = INTEGER(Rlind);
    int* crind = INTEGER(Rrind);
    int maxInd = 0;
    for(int i = 0; i < n; i++){ maxInd = max(maxInd, crind[i]); }
    baseCH.resize(maxInd+2);
    baseS.resize(maxInd+2);
    baseP.resize(maxInd+1);

	double denomMax = maxInd + 1.0;
	double startProb = 1.0 / denomMax;
	double tot = 0;
	for(int i = 0; i <= maxInd; i++){ 
		baseP[i] = startProb; 
		tot += startProb;
	}
		
	
	p2s();
	s2ch();
	
	int this_l, this_r;
	obs_inf.resize(n);
	node_inf.resize(maxInd+2);
	for(int i = 0; i < n; i++){
		this_l = clind[i];
		this_r = crind[i];
		obs_inf[i].l = this_l;
		obs_inf[i].r = this_r;
		node_inf[this_l].l.push_back(i);
		node_inf[this_r].r.push_back(i);
	}
	current_llk = R_NegInf;
}

void emicm::p2s(){
	int k = baseP.size();
	baseS.resize(k + 1);
	baseS[0] = 1.0;
	baseS[k] = 0.0;
	for(int i = 1; i < k; i++){ baseS[i] = baseS[i-1] - baseP[i-1]; }
}

void emicm::s2ch(){
	int k = baseS.size();
	baseCH.resize(k);
	baseCH[0] = R_NegInf;
	baseCH[k-1] = R_PosInf;
	for(int i = 1; i < (k-1); i++){ baseCH[i] = log(-log(baseS[i]));}
}

void emicm::ch2p(){
	int k = baseCH.size();
	baseS[0] = 1.0;
	baseS[k-1] = 0.0;
	for(int i = 1; i < (k-1); i++){ baseS[i] = exp(-exp(baseCH[i]));}
	for(int i = 1; i < k; i++){ baseP[i-1] = baseS[i-1] - baseS[i];}
}



