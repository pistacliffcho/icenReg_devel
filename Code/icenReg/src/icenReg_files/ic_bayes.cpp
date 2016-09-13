#include "ic_bayes.h"


//		POSTERIOR CALCULATION FUNCTIONS
double IC_bayes::computePriorLogDens(Eigen::VectorXd &propVec){
	Rcpp::NumericVector rPropVec = eigen2RVec(propVec);
	Rcpp::NumericVector r_ans = priorFxn(rPropVec);
	double ans = r_ans[0];
	return(ans);
}

double IC_bayes::computeLLK(Eigen::VectorXd &propVec){
	int propSize = propVec.size();
	int baselineSize = baseIC->b_pars.size();
	int regSize = baseIC->betas.size();
	if(propSize != (baselineSize + regSize) ){
		Rprintf("Error: propSize != regSize\n");
		return(0.0);
	}
	for(int i = 0; i < baselineSize; i++){
		baseIC->b_pars[i] = propVec[i];
	}
	for(int i = 0; i < regSize; i++){
		baseIC->betas[i] = propVec[i + baselineSize];
	}
	double ans = baseIC->calcLike_all();
	return(ans);
}

double IC_bayes::computePosteriorLogDens(Eigen::VectorXd &propVec){
	double ans = computePriorLogDens(propVec);
	ans += computeLLK(propVec);
	return(ans);
}



//          SAMPLING FUNCTIONs
void MHBlockUpdater::proposeNewParameters(){
	proposalParameters.resize(totParams);
	for(int i = 0; i < totParams; i++){
		proposalParameters[i] = R::rnorm(0.0, 1.0);
	}	
	proposalParameters = cholDecomp * proposalParameters.transpose() + currentParameters;
	proposeLogDens = logPostDens(proposalParameters, posteriorCalculator);
}

void MHBlockUpdater::updateCholosky(){
	Rprintf("Not currently working\n");
}

void MHBlockUpdater::acceptOrReject(){
	if(proposeLogDens >= currentLogDens){
		currentLogDens = proposeLogDens;
		currentParameters = proposalParameters;
	}
	else{
		double acceptProb = exp(proposeLogDens - currentLogDens);
		double randVal = R::runif(0.0, 1.0);
		if(randVal < acceptProb){
			currentLogDens = proposeLogDens;
			currentParameters = proposalParameters;
		}
	}
}

void MHBlockUpdater::mcmc(){
	
	if(logPostDens == NULL){
		throw std::range_error("logPostDens pointer not initialized in MHBlockUpdater.\n");
	}
	if(posteriorCalculator == NULL){
		throw std::range_error("posteriorCalculator not initialized in MHBlockUpdater.\n");
	}
	
	currentLogDens = R_NegInf;
	proposeNewParameters();
	acceptOrReject();

	
	savedValues.resize(numSaved, totParams);
	savedLPD.resize(numSaved);
	
	cout << "Starting parameters: \n" << currentParameters << 
			"\nStarting Chol: \n" << cholDecomp << "\n";
	
	
	int saveCount = 0;
	for(int i = 0; i < samples; i++){
		proposeNewParameters();
		acceptOrReject();
		if( (i % thin) == 0){
			savedValues.row(saveCount) = currentParameters;
			savedLPD[i]                = currentLogDens;
			saveCount++;
		}
		if( (i % iterationsPerUpdate) && updateChol ){
			updateCholosky();
		}
	}
}

double logIC_bayesPostDens(Eigen::VectorXd &propVec, void* void_icBayesPtr){
	IC_bayes* icBayesPtr = static_cast<IC_bayes*>(void_icBayesPtr);
	double ans = icBayesPtr->computePosteriorLogDens(propVec);
	return(ans);
}





//		CONSTRUCTOR AND DECONSTRUCTOR FUNCTIONS
MHBlockUpdater::MHBlockUpdater(Eigen::VectorXd &initValues, Eigen::MatrixXd &initChol,
				  int samples, int thin, int iterationsPerUpdate,
				  bool updateChol){
	currentParameters = initValues;
	cholDecomp = initChol;
	this->samples = samples;
	this->thin = thin;
	this->iterationsPerUpdate = iterationsPerUpdate;
	this->updateChol = updateChol;		
	double numSaved_double =((double) samples) / ((double) thin);
	numSaved = floor(numSaved_double);
	totParams = currentParameters.size();	
	logPostDens = NULL;	  
	posteriorCalculator = NULL;
}

IC_bayes::IC_bayes(Rcpp::List R_bayesList, Rcpp::Function R_priorFxn,
				   Rcpp::List R_ic_parList)
                   : priorFxn(R_priorFxn){
                   
    Rcpp::IntegerVector R_linkType = R_ic_parList["linkType"];
    int cLinkType = INTEGER(R_linkType)[0];
    if(cLinkType == 1 || cLinkType == 2){ 
	    baseIC = new IC_parOpt(R_ic_parList);
    }
    else if(cLinkType == 3){
    	baseIC = new IC_parOpt_aft(R_ic_parList);
    }
    else{
    	Rprintf("Warning: invalid link type! Setting to aft\n");
    	baseIC = new IC_parOpt_aft(R_ic_parList);
    }          
    int n_reg_pars = baseIC->betas.size();
    int n_b_pars   = baseIC->b_pars.size();
     
    Eigen::MatrixXd eChol;
    
    Rcpp::LogicalVector R_useMLE_start = R_bayesList["useMLE_start"];
    Rcpp::IntegerVector R_samples      = R_bayesList["samples"];
    Rcpp::IntegerVector R_thin         = R_bayesList["thin"];
    Rcpp::IntegerVector R_it_per_up    = R_bayesList["iterationsPerUpdate"];
    Rcpp::LogicalVector R_updateChol   = R_bayesList["updateChol"];

    baseIC->successfulBuild = true;
    if(Rf_isNull(R_useMLE_start)) baseIC->successfulBuild = false;
    if(Rf_isNull(R_samples))      baseIC->successfulBuild = false;
    if(Rf_isNull(R_thin))         baseIC->successfulBuild = false;
    if(Rf_isNull(R_it_per_up))    baseIC->successfulBuild = false;
    if(Rf_isNull(R_updateChol))   baseIC->successfulBuild = false;
    
    if(!baseIC->successfulBuild){ return; }
    
    bool useMLE_start        = LOGICAL(R_useMLE_start)[0] == TRUE;
    int samples              = INTEGER(R_samples)[0];
    int thin                 = INTEGER(R_thin)[0];
    int iterationsPerUpdate  = INTEGER(R_it_per_up)[0];
    bool updateChol          = LOGICAL(R_updateChol)[0] == TRUE;
    
    if(useMLE_start){
        baseIC->optimize();
            	 	
        // This is a bit backwards, but since the tools have
        // already been built to fill the full Hessian into 
        // an R-matrix
            	 	
        Rcpp::NumericVector R_score(n_reg_pars + n_b_pars);
        Rcpp::NumericMatrix R_Hessian(n_reg_pars + n_b_pars);
		baseIC->fillFullHessianAndScore(R_Hessian, R_score);
		Eigen::MatrixXd eHess;
		copyRmatrix_intoEigen(R_Hessian, eHess);
		eHess = -eHess.inverse();   
		eChol =  eHess.llt().matrixL();       	 	
    }
    else{
    	Rprintf("NOTE: ADAPTIVE BLOCK UPDATING IS NOT PROPERLY IMPLEMENTED YET!!!\n");
    	eChol.resize(n_reg_pars + n_b_pars, n_reg_pars + n_b_pars);
    	eChol *= 0.0;
    	for(int i = 0; i < (n_reg_pars + n_b_pars); i++){
    		eChol(i, i) = 1.0;
    	}
	}
	Eigen::VectorXd initPars(n_b_pars + n_reg_pars);
	for(int i = 0; i < n_b_pars; i++){
		initPars[i] = baseIC->b_pars[i];
	}
	for(int i = 0; i < n_reg_pars; i++){
		initPars[i + n_b_pars] = baseIC->betas[i];
	}
	mcmcInfo = new MHBlockUpdater(initPars, eChol,
				  samples, thin, iterationsPerUpdate,
				  updateChol);
	
	mcmcInfo->logPostDens         = logIC_bayesPostDens;
	mcmcInfo->posteriorCalculator = this;
}

IC_bayes::~IC_bayes(){
	delete mcmcInfo;
	delete baseIC;
}



//      R Interface 
Rcpp::List R_ic_bayes(Rcpp::List R_bayesList, Rcpp::Function priorFxn, 
					  Rcpp::List R_ic_parList){
	

    IC_bayes bayes(R_bayesList, priorFxn, R_ic_parList);          
	
	if(!bayes.baseIC->successfulBuild){
		Rprintf("Unsuccessful build of C++ IC_bayes object!\n");
		Rcpp::List ans;
		return(ans);
	}
	
	bayes.mcmcInfo->mcmc();	
	
	Rcpp::List ans;
	ans["samples"] = eigen2RMat(bayes.mcmcInfo->savedValues);
	ans["logPosteriorDensity"] = eigen2RVec(bayes.mcmcInfo->savedLPD);
	return(ans);
}
