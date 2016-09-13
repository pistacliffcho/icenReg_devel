class MHBlockUpdater{
	public : 
	double (*logPostDens)(Eigen::VectorXd &parameters, void* otherInfo);
	// pointer to a function that calculates log posterior density 
	// for a given set of parameters
	
	int samples, thin, iterationsPerUpdate, totParams, numSaved;
	void *posteriorCalculator;
	bool updateChol;
	double currentLogDens, proposeLogDens;
	
	Eigen::VectorXd currentParameters, proposalParameters, savedLPD;
	Eigen::MatrixXd savedValues;
	Eigen::MatrixXd cholDecomp;	
	
	void proposeNewParameters();
	void acceptOrReject(); 
	void updateCholosky();
	void mcmc();
	
	MHBlockUpdater(Eigen::VectorXd &initValues, Eigen::MatrixXd &initChol,
				  int samples, int thin, int iterationsPerUpdate,
				  bool updateChol);	
};


class IC_bayes {
	public :
	
	
	IC_parOpt* baseIC;
	Rcpp::Function priorFxn;
	double computePriorLogDens(Eigen::VectorXd &propVec);
	double computeLLK(Eigen::VectorXd &propVec);
	double computePosteriorLogDens(Eigen::VectorXd &propVec);
	MHBlockUpdater* mcmcInfo;
	
	IC_bayes(Rcpp::List R_bayesList, Rcpp::Function R_priorFxn,
			  Rcpp::List R_ic_parList);
	
    ~IC_bayes();
};

double logIC_bayesPostDens(Eigen::VectorXd &propVec, void* void_icBayesPtr);

