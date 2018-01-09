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
	baseIC->update_etas();
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
	proposalParameters = cholDecomp * proposalParameters + currentParameters;
	proposeLogDens = logPostDens(proposalParameters, posteriorCalculator);
}

void MHBlockUpdater::updateCholesky(Eigen::MatrixXd valMat){
	int nCols = valMat.cols();
	timesAdapted++;
	double acceptRate = timesAccepted / timesRan;
	gamma1 = 1.0/pow(timesAdapted + 3.0, 0.8);
    double gamma2 = 10.0 * gamma1;
    double adaptFactor = exp(gamma2 * (acceptRate - optimalAR) );
    cholScale = cholScale * adaptFactor;
		
	double colMean;
	double denom = valMat.rows();
	int nRows = valMat.rows();
	for(int j = 0; j < nCols; j++){
		colMean = 0;
		for(int i = 0; i < nRows; i++){ colMean += valMat(i, j); }
		colMean = colMean / denom;
		for(int i = 0; i < nRows; i++){
			valMat(i,j) -= colMean;
		}	
	}
	Eigen::MatrixXd intermMat = xtx(valMat, 0, nRows - 1) / (denom - 1.0);
	
	propCov = propCov + gamma1 * (intermMat - propCov);
	cholDecomp = propCov.llt().matrixL();
	
	timesRan = 0;
	timesAccepted = 0;
}

void MHBlockUpdater::acceptOrReject(){
	timesRan++;
	if(ISNAN(proposeLogDens)){return;}
	if(proposeLogDens >= currentLogDens){
		currentLogDens = proposeLogDens;
		currentParameters = proposalParameters;
		timesAccepted++;
	}
	else{
		double acceptProb = exp(proposeLogDens - currentLogDens);
		double randVal = R::runif(0.0, 1.0);
		if(randVal < acceptProb){
			currentLogDens = proposeLogDens;
			currentParameters = proposalParameters;
			timesAccepted++;
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
	
	timesRan      = 0;
	timesAccepted = 0;
	timesAdapted  = 0;
	
	currentLogDens = logPostDens(currentParameters, posteriorCalculator);
	proposeNewParameters();
	acceptOrReject();
	
	int numberParms = currentParameters.size();
	burnInMat.resize(burnIn, numberParms );
	burnInMat *= 0.0;
	for(int j = 0; j < burnIn/2; j++){
		proposeNewParameters();
		acceptOrReject();		
		burnInMat.row(j) = currentParameters;	
	}
	
	for(int j = 0; j < burnIn/2; j++){
	  proposeNewParameters();
	  acceptOrReject();		
	  burnInMat.row(j) = currentParameters;	
  	int start_row, stop_row;
  	if(updateChol & ( (j + 1) % samplesPerUpdate == 0) ){
  	  start_row = j - samplesPerUpdate + 1;
  	  stop_row  = j - 1; 
  	  Eigen::MatrixXd recentVals = copyRows(burnInMat, 
                                            start_row, 
                                            stop_row); 
	    updateCholesky(recentVals);  
	  }
	}
	
	
		
	savedValues.resize(samples, totParams);
	savedLPD.resize(samples);
	
	for(int i = 0; i < samples; i++){
		for(int j = 0; j < thin; j++){
			proposeNewParameters();
			acceptOrReject();
		}
		savedValues.row(i) = currentParameters;
		savedLPD[i]        = currentLogDens;
	}
}

double logIC_bayesPostDens(Eigen::VectorXd &propVec, void* void_icBayesPtr){
	IC_bayes* icBayesPtr = static_cast<IC_bayes*>(void_icBayesPtr);
	double ans = icBayesPtr->computePosteriorLogDens(propVec);
	return(ans);
}





//		CONSTRUCTOR AND DECONSTRUCTOR FUNCTIONS
MHBlockUpdater::MHBlockUpdater(Eigen::VectorXd &initValues, Eigen::MatrixXd &initCov,
				  int samples, int thin, int samplesPerUpdate,
				  bool updateChol, int burnIn, double cholScale, 
				  double acceptRate){
	currentParameters         = initValues;
	propCov                   = initCov;
	cholDecomp                = propCov.llt().matrixL();
	this->samples             = samples;
	this->thin                = thin;
	this->samplesPerUpdate    = samplesPerUpdate;
	this->updateChol          = updateChol;		
	this->burnIn              = burnIn;
	double numSaved_double    = ((double) samples) / ((double) thin);
	this->cholScale           = cholScale;
	numSaved                  = ceil(numSaved_double);
	totParams                 = currentParameters.size();	
	logPostDens               = NULL;	  
	posteriorCalculator       = NULL;
	timesRan                  = 0.0;
	timesAccepted             = 0.0;
	gamma1                    = 0.0;
	optimalAR                 = acceptRate;
	timesAdapted              = 0.0;
		
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
    int n_reg_pars  = baseIC->betas.size();
    int n_b_pars    = baseIC->b_pars.size();
    int totalParams = n_reg_pars + n_b_pars;

//    Eigen::MatrixXd eChol;
    
    Rcpp::LogicalVector R_useMLE_start = R_bayesList["useMLE_start"];
    Rcpp::IntegerVector R_samples      = R_bayesList["samples"];
    Rcpp::IntegerVector R_thin         = R_bayesList["thin"];
    Rcpp::IntegerVector R_sample_per_up    = R_bayesList["samplesPerUpdate"];
    Rcpp::LogicalVector R_updateChol   = R_bayesList["updateChol"];
    Rcpp::NumericVector R_cholScale    = R_bayesList["initSD"];
    Rcpp::IntegerVector R_burnIn       = R_bayesList["burnIn"];
    Rcpp::NumericVector R_acceptRate   = R_bayesList["acceptRate"];
    
    bool useMLE_start        = LOGICAL(R_useMLE_start)[0] == TRUE;
    int samples              = INTEGER(R_samples)[0];
    int thin                 = INTEGER(R_thin)[0];
    int samplesPerUpdate     = INTEGER(R_sample_per_up)[0];
    double acceptRate        = REAL(R_acceptRate)[0];
    bool updateChol          = LOGICAL(R_updateChol)[0] == TRUE;
    double cholScale         = R_cholScale[0];
    int burnIn               = R_burnIn[0];

	Eigen::MatrixXd propCov;

    if(useMLE_start){
        baseIC->optimize();
            	 	
        // This is a bit backwards, but since the tools have
        // already been built to fill the full Hessian into 
        // an R-matrix
            	 	
        Rcpp::NumericVector R_score(n_reg_pars + n_b_pars);
        Rcpp::NumericMatrix R_Hessian(n_reg_pars + n_b_pars);
		baseIC->fillFullHessianAndScore(R_Hessian, R_score);
		copyRmatrix_intoEigen(R_Hessian, propCov);
		Eigen::MatrixXd intermMat = -propCov.inverse();   
		propCov = intermMat;
    }
    else{
    	propCov.resize(totalParams, totalParams);
//    	propCov *= 0.0;
        for(int i = 0; i < (totalParams); i++){
          for(int j = 0; j < (totalParams); j++){
            propCov(i,j) = 0.;
          }  
        }
    	for(int i = 0; i < (totalParams); i++){
    		propCov(i, i) = cholScale;
    	}	    	
	}
		
	Eigen::VectorXd initPars(n_b_pars + n_reg_pars);
	if(useMLE_start){
		for(int i = 0; i < n_b_pars; i++){	
			initPars[i] = baseIC->b_pars[i];
		}
		for(int i = 0; i < n_reg_pars; i++){
			initPars[i + n_b_pars] = baseIC->betas[i];
		}
	}
	else{
		for(int i = 0; i < totalParams; i++){ initPars[i] = 0; }
	}
	mcmcInfo = new MHBlockUpdater(initPars, propCov,
				  samples, thin, samplesPerUpdate,
				  updateChol, burnIn, cholScale, acceptRate);
	
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
	ans["finalChol"] = eigen2RMat(bayes.mcmcInfo->propCov);
	return(ans);
}
