class MHBlockUpdater{
	double (*logPostDens)(Eigen::VectorXd &parameters, void* otherInfo);
	// pointer to a function that calculates log posterior density 
	// for a given set of parameters
	
	int samples, thin, iter, iterationsPerUpdate, totParams;
	void *posteriorCalculator;
	bool updateChol;
	double currentLogDens, proposeLogDens;
	
	Eigen::VectorXd currentParameters, proposalParameters;
	Eigen::MatrixXd savedValues;
	Eigen::MatrixXd cholDecomp;
	
	void proposeNewParameters();
	void acceptOrReject(); // Makes the decision to keep
	void updateCholosky();
	void mcmc(int samples);	
};

class IC_bayes : public IC_parOpt{
	Rcpp::Function priorFxn;
	double computePriorLogDens(Eigen::VectorXd &propVec);
	double computeLLK(Eigen::VectorXd &propVec);
	double computePosteriorLogDens(Eigen::VectorXd &propVec);
	MHBlockUpdater mcmcInfo;
	
	IC_bayes(SEXP R_prior, SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w);
};

double logIC_bayesPostDens(Eigen::VectorXd &propVec, void* icBayesPtr);

