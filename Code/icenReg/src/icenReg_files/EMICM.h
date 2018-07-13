class emicm{
public:
    void update_p_ob(int i, bool useS);  
    void update_p_ob(bool useS); 
    double llk(bool useS);           
    double current_llk;
    double tot_w;
    
	Eigen::VectorXd baseP;
	Eigen::VectorXd baseS;
	Eigen::VectorXd baseCH;
	Eigen::VectorXd innerCH;
	Eigen::VectorXd pobs;
  double* w;
    
	void p2s();
	void s2ch();
	void ch2p();

	Eigen::VectorXd em_m; 
	Eigen::VectorXd ch_d1;
	Eigen::VectorXd ch_d2;    
	
	Eigen::VectorXd ch_d1_con;
	Eigen::VectorXd ch_d2_con1;
	Eigen::VectorXd ch_d2_con2;
	Eigen::VectorXd prop_delta;
	
    vector<obInf> obs_inf;
    vector<node_info> node_inf;
    
    emicm(SEXP Rlind, SEXP Rrind, SEXP R_w);
    void em_step(int iters);
    void icm_step();
    
    void calc_m_for_em();
    void calc_icm_ders();
    
    double run(double tol, int maxIter, int emSteps);
    
   	int iter; 
};


void addIcmProp(Eigen::VectorXd &bch, Eigen::VectorXd &prop);
double icmFirstDerv(double prt1, double pob);
double icmSecondDerv(double prt2, double prt3, double pob, bool isLeft);
extern "C" {
SEXP EMICM(SEXP Rlind, SEXP Rrind, SEXP iters, SEXP R_w);
}