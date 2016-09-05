class MHBlockUpdater{
	public : 
	double (*logPostDens)(Eigen::VectorXd &parameters, void* otherInfo);
	// pointer to a function that calculates log posterior density 
	// for a given set of parameters
	
	int samples, thin, iterationsPerUpdate, totParams, numSaved;
	void *posteriorCalculator;
	bool updateChol;
	double currentLogDens, proposeLogDens;
	
	Eigen::VectorXd currentParameters, proposalParameters;
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


class IC_bayes : public IC_parOpt{
	public :
	Rcpp::Function priorFxn;
	double computePriorLogDens(Eigen::VectorXd &propVec);
	double computeLLK(Eigen::VectorXd &propVec);
	double computePosteriorLogDens(Eigen::VectorXd &propVec);
	MHBlockUpdater* mcmcInfo;
	
	IC_bayes(Rcpp::Function R_prior, SEXP useMLE_start, SEXP updateChol,
			  SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w);
    ~IC_bayes();
};

double logIC_bayesPostDens(Eigen::VectorXd &propVec, void* void_icBayesPtr);

