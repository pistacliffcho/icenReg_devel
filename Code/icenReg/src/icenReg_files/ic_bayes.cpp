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
