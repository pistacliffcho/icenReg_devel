class emicm{
public:
    void update_p_ob(int i, bool useS);  
    void update_p_ob(bool useS); 
    double llk(bool useS);           
    double current_llk;
    
	Eigen::VectorXd baseP;
	Eigen::VectorXd baseS;
	Eigen::VectorXd baseCH;
	Eigen::VectorXd pobs;

	void p2s();
	void s2ch();
	void ch2p();

	Eigen::VectorXd em_m; 
	Eigen::VectorXd ch_d1;
	Eigen::VectorXd ch_d2;    
	
    vector<obInf> obs_inf;
    vector<node_info> node_inf;
    
    emicm(SEXP Rlind, SEXP Rrind);
    void em_step();
    void icm_step();
    
    void calc_m_for_em();
};

void calc_m_for_em(){
	double this_m = 0;
	node_info* thisNode;
	int k = node_inf.size();
	double n = obs_inf.size();
	for(int i = 0; i < k; i++){
		thisNode = &node_inf[i];
		for(unsigned int j = 0; j < thisNode->l.size(); j++){
			this_m += 1.0 / pobs[ thisNode->l[j] ]
		} 
		for(unsigned int j = 0; j < thisNode->r.size(); j++){
			this_m -= 1.0 / pobs[ thisNode->r[j] ]
		}
		em_m[i] = this_m / n; 
	}
}

void em_step(){
	p2s();
	update_p_ob(true);
	calc_m_for_em();
	int k = baseP.size();
	calc_m_for_em();
	for(int j = 0; j < k; j++){ baseP[i] *= em_m[i];}
	p2s();
	update_p_ob(true);
}

void emicm::update_p_ob(int i, bool useS){
	if(useS){ pobs[i] = baseS[ obs_inf[i].l ] - baseS[ obs_inf[i].r + 1]; }
	else{
		double chl = baseCH[ obs_inf[i].l ] 
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
	return(current_llk);
}

emicm::emicm(SEXP Rlind, SEXP Rrind){
    int n = LENGTH(Rlind);
    if(n != LENGTH(Rrind)){Rprintf("length of Rlind and Rrind not equal\n"); return;
    
    pobs.resize(n);
    
    int* clind = INTEGER(Rlind);
    int* crind = INTEGER(Rrind);
    int maxInd = 0;
    for(int i = 0; i < n; i++){ maxInd = max(maxInd, crind[i]); }
    baseCH.resize(maxInd+2);
    baseS.resize(maxInd+2);
    baseP.resize(maxInd+1);

	double denomMax = maxInd + 1.0;
	double startProb = 1 / denomMax;
	for(int i = 0; i <= maxInd; i++){ baseP[i] = startProb; }
	p2s();
	s2ch();
	
	int this_l, this_r;
	obs_inf.resize(n);
	node_inf.resize(maxInd+2);
	vector<int> minActPoints;
	for(int i 0; i < n; i++){
		this_l = clind[i];
		this_r = crind[i];
		addIfNeeded(minActPoints, this_l, this_r, maxInd);
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
	baseCH[0] = 0;
	baseCH[k] = R_PosInf;
	for(int i = 1; i < (k-1); i++){ baseCH[i] = log(-log(baseS[i]));}
}

void emicm::ch2p(){
	int k = baseCH.size();
	baseS[0] = 1.0;
	baseS[k-1] = 0.0;
	for(int i = 1; i < (k-1); i++){ baseS[i] = exp(-exp(baseCH[i]));}
	for(int i = 1; i < k; i++){ baseP[i-1] = baseS[i-1] - baseS[i];}
}
