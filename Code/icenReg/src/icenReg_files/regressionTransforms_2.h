// Changes to class:
//   Always uses Rcpp types
//   Should be fed baseline parameters
//   Handles AFT model internally rather than relying on native R code to handle
class condProbCal_2{
public:
    double (*getBase_s)(double q, vector<double> &bli);
    double (*getBase_q)(double s, vector<double> &bli);
    double (*baseSurv_2_condSurv)(double s, double nu);
    double (*condSurv_2_baseSurv)(double s, double nu);
//    double (*transform_p)(double q, double nu);
    double calc_p(double q, double nu, vector<double> &bli);
    double calc_q(double p, double nu, vector<double> &bli);
    vector<double> preppedParams;
	condProbCal_2(Rcpp::CharacterVector regType, 
				  Rcpp::CharacterVector baseType);
	bool isOK, isAFT;
};

std::vector<double> getRow(int row, Rcpp::NumericMatrix rMat){
	int nCol = rMat.ncol();
	int nRow = rMat.nrow();
	std::vector<double> ans(nCol);
	if(row < nRow){
		for(int j = 0; j < nCol; j++){ ans[j] = rMat(row, j); }
	}
	return(ans);
}

// SURVIVAL TRANSFORMS
double baseSurv_2_condSurv_ph(double s, double nu){
	if(s == 0 || s == 1) return(s);
    double ans = pow(s, nu);
    return( ans );
}

double baseSurv_2_condSurv_po(double s, double nu){
	if(s == 0 || s == 1) return(s);
    double ans;
    double prod = s * nu;
    ans = prod/(prod - s + 1);
    return(ans);
}

double condSurv_2_baseSurv_ph(double s, double nu){
	if(s == 0.0 || s == 1.0) return(s);
    double log_s = log(s);
    double log_trans = log_s / nu;
    return( exp(log_trans) );
}
double condSurv_2_baseSurv_po(double s, double nu){
	if(s == 1.0 || s == 0.0) return(s);
    return( s * (1.0/nu) / (s * 1.0/nu - s + 1.0) );
}


// BASELINE MODELS
double getExpSurv(double q, vector<double> &bli){
    return( R::pexp(q, exp(bli[0]), 0, 0) );
}
double getWeibSurv(double q, vector<double> &bli){
    return( R::pweibull(q, exp(bli[0]), exp(bli[1]), 0, 0));
}
double getLogNormSurv(double q, vector<double> &bli){
    return( R::pnorm(log(q), bli[0], exp(bli[1]), 0, 0) );
}
double getGammaSurv(double q, vector<double> &bli){
    return( R::pgamma(q, exp(bli[0]), exp(bli[1]), 0,0 ) );
}
double getLgLgsticSurv(double q, vector<double> &bli){
    return( 1.0 - ic_ploglogistic(q, exp(bli[0]), exp(bli[1])));
}
double getGenGammaSurv(double q, vector<double> &bli){
    return(1.0 - ic_pgeneralgamma(q, bli[0], exp(bli[1]), bli[2]));
}




double getExpQ(double s, vector<double> &bli){ return( R::qexp(1.0 - s, exp(bli[0]), 1, 0) ); }

double getWeibQ(double s, vector<double> &bli){
    return( R::qweibull(1.0-s, exp(bli[0]), exp(bli[1]), 1, 0));
}
double getLogNormQ(double s, vector<double> &bli){
    return( exp( R::qnorm(1.0-s, bli[0], exp(bli[1]), 1, 0) ) );
}
double getGammaQ(double s, vector<double> &bli){
    return( R::qgamma(1.0-s, exp(bli[0]), exp(bli[1]), 1,0 ) );
}
double getLgLgsticQ(double s, vector<double> &bli){
    return(ic_qloglogistic(1.0-s, exp(bli[0]), exp(bli[1]) ) );
}
double getGenGammaQ(double s, vector<double> &bli){
    return(ic_qgeneralgamma(1.0-s, bli[0], exp(bli[1]), bli[2] ) );
}




