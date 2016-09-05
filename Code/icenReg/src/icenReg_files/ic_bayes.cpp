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

IC_bayes::IC_bayes(Rcpp::Function R_prior, SEXP useMLE_start, SEXP updateChol,
			  SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w) 
              : IC_parOpt(R_s_t, R_d_t, R_covars,
              R_uncenInd, R_gicInd, R_lInd, R_rInd,
              R_parType, R_linkType, R_w),
              priorFxn(R_prior)
              {
}

IC_bayes::~IC_bayes(){
	delete mcmcInfo;
}