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
	int baselineSize = b_pars.size();
	int regSize = betas.size();
	if(propSize != (baselineSize + regSize) ){
		Rprintf("Error: propSize != regSize\n");
		return(0.0);
	}
	for(int i = 0; i < baselineSize; i++){
		b_pars[i] = propVec[i];
	}
	for(int i = 0; i < regSize; i++){
		betas[i] = propVec[i + baselineSize];
	}
	double ans = calcLike_all();
	return(ans);
}

double IC_bayes::computePosteriorLogDens(Eigen::VectorXd &propVec){
	double ans = computePriorLogDens(propVec);
	ans += computeLLK(propVec);
	return(ans);
}



//          SAMPLING FUNCTIONs
void MHBlockUpdater::proposeNewParameters(){
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
		double acceptProb = exp(currentLogDens - proposeLogDens);
		double randVal = R::runif(0.0, 1.0);
		if(randVal < acceptProb){
			currentLogDens = proposeLogDens;
			currentParameters = proposalParameters;
		}
	}
}

void MHBlockUpdater::mcmc(){
	savedValues.resize(numSaved, totParams);
	int saveCount = 0;
	for(int i = 0; i < samples; i++){
		proposeNewParameters();
		acceptOrReject();
		if( (i % thin) == 0){
			savedValues.row(saveCount) = currentParameters;
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
}

IC_bayes::IC_bayes(Rcpp::List R_bayesList, Rcpp::Function R_priorFxn,
			  SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w) 
	:IC_parOpt(R_s_t, R_d_t, R_covars,
               R_uncenInd, R_gicInd, R_lInd, R_rInd,
               R_parType, R_linkType, R_w),
     priorFxn(R_priorFxn){
    int n_reg_pars = betas.size();
    int n_b_pars = b_pars.size(); 
    Eigen::MatrixXd eChol;
    
    Rcpp::LogicalVector R_useMLE_start = R_bayesList["useMLE_start"];
    Rcpp::IntegerVector R_samples      = R_bayesList["samples"];
    Rcpp::IntegerVector R_thin         = R_bayesList["thin"];
    Rcpp::IntegerVector R_it_per_up    = R_bayesList["iterationsPerUpdate"];
    Rcpp::LogicalVector R_updateChol   = R_bayesList["updateChol"];
    
    bool useMLE_start        = LOGICAL(R_useMLE_start)[0] == TRUE;
    int samples              = INTEGER(R_samples)[0];
    int thin                 = INTEGER(R_thin)[0];
    int iterationsPerUpdate  = INTEGER(R_it_per_up)[0];
    bool updateChol          = LOGICAL(R_updateChol)[0] == TRUE;
    
    if(useMLE_start){
        optimize();
            	 	
        // This is a bit backwards, but since the tools have
        // already been built to fill the full Hessian into 
        // an R-matrix
            	 	
        Rcpp::NumericVector R_score(n_reg_pars + n_b_pars);
        Rcpp::NumericMatrix R_Hessian(n_reg_pars + n_b_pars);
		fillFullHessianAndScore(R_Hessian, R_score);
		Eigen::MatrixXd eHess;
		copyRmatrix_intoEigen(R_Hessian, eHess);
		eChol = eHess.llt().matrixL();            	 	
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
		initPars[i] = b_pars[i];
	}
	for(int i = 0; i < n_reg_pars; i++){
		initPars[i + n_b_pars] = betas[i];
	}
	mcmcInfo = new MHBlockUpdater(initPars, eChol,
				  samples, thin, iterationsPerUpdate,
				  updateChol);
	
}

IC_bayes::~IC_bayes(){
	delete mcmcInfo;
}