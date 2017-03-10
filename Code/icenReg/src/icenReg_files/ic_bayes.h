class MHBlockUpdater{
	public : 
	double (*logPostDens)(Eigen::VectorXd &parameters, void* otherInfo);
	// pointer to a function that calculates log posterior density 
	// for a given set of parameters
	
	int samples, thin, samplesPerUpdate, totParams, numSaved,
		burnIn;
	void *posteriorCalculator;
	bool updateChol;
	double currentLogDens, proposeLogDens, cholScale,
			timesRan, timesAccepted, timesAdapted, gamma1, optimalAR;
	
	Eigen::VectorXd currentParameters, proposalParameters, savedLPD;
	Eigen::MatrixXd savedValues;
	Eigen::MatrixXd propCov;
	Eigen::MatrixXd cholDecomp;	
	Eigen::MatrixXd burnInMat;

	
	
	void proposeNewParameters();
	void acceptOrReject(); 
	void updateCholesky(Eigen::MatrixXd valMat);
	void mcmc();
	
	MHBlockUpdater(Eigen::VectorXd &initValues, Eigen::MatrixXd &initCov,
				  int samples, int thin, int samplesPerUpdate,
				  bool updateChol, int burnIn, double cholScale, 
				  double acceptRate);	
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

