//
//  ic_par.h
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#ifndef ____ic_par__
#define ____ic_par__
class linkFun{
//Abstract class with necessary information about link function
public:
    virtual double con_s(double b_s, double nu) =0;
    //Computes conditional survival probability from baseline S and nu
    virtual double con_d(double b_d, double b_s, double nu) =0;
    //Computes condtional density from baseline density, baseline S and nu
    
/*    virtual double con_s_der(double b_s, double nu) = 0;
    //Computes derivative of conditional survival probability with respect to nu
    virtual double con_d_der(double b_d, double b_s, double nu) = 0;
    //Computes derviative of conditional density with respect to nu
    virtual double con_s_der2(double b_s, double nu) = 0;
    //Computes 2nd derivative of conditional survival probability with respect to nu
    virtual double con_d_der2(double b_d, double b_s, double nu) = 0;
    //Computes 2nd derviative of conditional density with respect to nu			*/
};

class propOdd : public linkFun{
// Proportional Odds link function
public:
    double con_s(double b_s, double nu) { return(nu * b_s / (b_s * (nu -1) + 1));}
    double con_d(double b_d, double b_s, double nu){
        double sqrt_denom = b_s * nu - b_s + 1;
        return( b_d * nu / (sqrt_denom * sqrt_denom));
    }
/*    double con_d_der(double b_d, double b_s, double nu){
        double ans = -b_d * (b_s * nu + b_s - 1)/pow(b_s * (nu - 1) + 1, 3.0);
        return(ans);
    }
    double con_s_der(double b_s, double nu){
        double ans = ((b_s - 1) * b_s) / pow( (b_s * (nu-1) + 1), 2.0);
        return(ans);
    }
    double con_s_der2(double b_s, double nu){
        double num = 2 * (b_s - 1) * b_s * b_s;
        double denom = pow(b_s * (nu - 1) + 1, 3.0);
        return(num/denom);
    }
    double con_d_der2(double b_d, double b_s, double nu){
        double num = 2 * b_s * b_s * ( b_s * (nu + 2) -2 );
        double denom = pow(b_s * (nu -1) + 1, 4.0);
        return(num/denom);
    }				*/
    virtual ~propOdd(){};
};

class propHaz : public linkFun{
// Proportional Hazards link function
public:
    double con_s(double b_s, double nu) { return(pow(b_s, nu));}
    double con_d(double b_d, double b_s, double nu){return( b_d * nu * pow(b_s, nu -1));}
/*    double con_d_der(double b_d, double b_s, double nu){
        double ans = b_d * pow(b_s, nu - 1);
        ans *= (nu * log(b_s) + 1);
        return(ans);
    }
    double con_s_der(double b_s, double nu){
        double ans = log(b_s) * pow(b_s, nu);
        return(ans);
    }
    double con_s_der2(double b_s, double nu){
        double logS = log(b_s);
        return(pow(b_s, nu) * logS * logS);
    }
    double con_d_der2(double b_d, double b_s, double nu){
        double logS = log(b_s);
        return(b_d * pow(b_s, nu - 1) * logS * (nu * logS + 2) );
    }				*/
    virtual ~propHaz(){};
};


class baseLineInfo{
//Abstract class with necessary information about baseline distirbution
public:
    virtual double base_d(double x, Eigen::VectorXd &par)=0;
    //computes baseline density
    virtual double base_s(double x, Eigen::VectorXd &par)=0;
    //computes baseline survival
    virtual void update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                                      Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals, Eigen::VectorXd &par)=0;
    //updates the baseline survival probabilities and densities required
};

class parBLInfo : public baseLineInfo{
//Abstract class with necessary info about parametric baseline distribution
    //Really just a placeholder: at some point I would like to add a spline-based baseline distribution
public:
    void update_baseline_vals(Eigen::VectorXd &s_t, Eigen::VectorXd &d_t,
                              Eigen::VectorXd &s_vals, Eigen::VectorXd &d_vals,
                              Eigen::VectorXd &par);
};

class gammaInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(R::dgamma(x, exp(par[0]), exp(par[1]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(R::pgamma(x, exp(par[0]), exp(par[1]), 0, 0));}
    virtual ~gammaInfo(){};
};

class weibullInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(R::dweibull(x, exp(par[0]), exp(par[1]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(R::pweibull(x, exp(par[0]), exp(par[1]), 0, 0));}
    virtual ~weibullInfo(){};
};

class lnormInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(ic_dlnorm(x, par[0], exp(par[1])));}
    double base_s(double x, Eigen::VectorXd &par){return(ic_plnorm(x, par[0], exp(par[1])));}
    virtual ~lnormInfo(){};
};

class expInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(R::dexp(x, exp(par[0]), 0));}
    double base_s(double x, Eigen::VectorXd &par){return(R::pexp(x, exp(par[0]), 0, 0));}
    virtual ~expInfo(){};
};

class loglogisticInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(ic_dloglogistic(x, exp(par[0]), exp(par[1])));}
    double base_s(double x, Eigen::VectorXd &par){return(1 - ic_ploglogistic(x, exp(par[0]), exp(par[1])));}
    virtual ~loglogisticInfo(){};
};

class genGammaInfo : public parBLInfo{
public:
    double base_d(double x, Eigen::VectorXd &par){return(ic_dgeneralgamma(x, par[0], exp(par[1]), exp(par[2]))) ;}
    double base_s(double x, Eigen::VectorXd &par){return(1.0 -ic_pgeneralgamma(x, par[0], exp(par[1]), exp(par[2]))) ;}
    virtual ~genGammaInfo(){};
};


struct dinf{
//structure for containing all the necessary info for uncensored data
    //values d, s and nu are the necessary indicators for looking up the density, survival and nu
public:
    int d;
    int s;
    int nu;
};

struct intInf{
// for general interval censoring
public:
    int l;
    int r;
    int nu;
};

struct rinf{
    //for right censored data
public:
    int l;
    int nu;
};

struct linf{
    //for left censored data
public:
    int r;
    int nu;
};

class IC_parOpt{
public:
    
    baseLineInfo* blInf;
    // Information about the baseline distribution
    
    linkFun* lnkFn;
    // Information about link function
    
    vector<double> w;                   //weights
    
    Eigen::VectorXd b_pars;                 //initialized
    //baseline distribution parameters
    Eigen::VectorXd d_b_pars;               //does not need initialization
    //derivatives of baseline parameters
    Eigen::MatrixXd d2_b_pars;              //does not need initialization
    //Hessian of baseline parameters
    Eigen::VectorXd betas;                  //initialized
    //regression parameters
    Eigen::VectorXd d_betas;                //does not need initialization
    // derivatives of betas
    Eigen::MatrixXd d2_betas;               //does not need initialization
    //Hessian of regression parameters
    Eigen::MatrixXd covars;                 //initialized
    //covariates
    Eigen::VectorXd eta;                    //initialized
    //vector of linear predictions
    Eigen::VectorXd expEta;                 //initialized
    //vector of exponential of linear predictors
    Eigen::VectorXd dobs_deta;
    //vector of derivatives of observation with respect to nu
    Eigen::VectorXd d2obs_d2eta;
    //vector of 2nd derviatives with respect to nu
    
    
    Eigen::VectorXd s_t;                    //initialized
    //vector of times associated with survival probabilities to calculate
    Eigen::VectorXd d_t;                    //initialized
    // vector of times associated with densities to caculate
    Eigen::VectorXd s_v;                    //Initialized at first call optObj.calcLike_all(), happens right outside constructor for IC_parOpt
    //vector of survival probabilities calculated at s_t
    Eigen::VectorXd d_v;
    // vector of densities calculated at d_t
    
    vector<dinf> uc;
    // vector of indicator information about uncensored observations
    // i.e. i's element is a pair of indicators for s_t and d_t values for subject i
    vector<intInf> gic;
    // same as above, but for general interval censored, i.e. l ind and r ind
    vector<linf> lc;
    // same as above, but for left censoring (i.e. indicator for right side of interval)
    vector<rinf> rc;
    // same as above, but for right censoring (i.e. indicator for left side of interval)

    virtual double calcLike_baseReady();
    //calculates the likelihood. Assumes baseline probs are ready
    virtual void calculate_baseline_probs(){
        blInf->update_baseline_vals(s_t,d_t,s_v, d_v,b_pars);};
    double calcLike_all(){
        calculate_baseline_probs();
        double ans = calcLike_baseReady();
        return(ans);
    };
    
    void NR_baseline_pars();
    void NR_reg_pars();

    void update_etas();
    virtual void update_dobs_detas();
    
    void calc_baseline_dervs();
    void numericCovar_dervs();
    void partAnalyticCovar_dervs();
    
    void fillFullHessianAndScore(SEXP r_mat, SEXP score);
    //for filling out the full hessian and score at the MLE

    double h, lk_new, lk_old;
	int iter, parType, linkType;
	bool successfulBuild;
	
	void optimize();

	// Initialization code to be shared by the two constructors
	void init(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w);
              
    Rcpp::List exportAns();
    
    IC_parOpt(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
              SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
              SEXP R_parType, SEXP R_linkType, SEXP R_w);
    IC_parOpt(Rcpp::List R_list);
    IC_parOpt(){}
    virtual ~IC_parOpt();
};

extern "C" {
    SEXP dGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q);
    SEXP pGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q);
    SEXP qGeneralGamma(SEXP R_x, SEXP R_mu, SEXP R_s, SEXP R_Q);
}


#endif /* defined(____ic_par__) */
