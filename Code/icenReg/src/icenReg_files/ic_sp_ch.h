//
//  ic_sp_ch.h
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#ifndef ____ic_sp_ch__
#define ____ic_sp_ch__
/*#include "../Eigen_local/Dense"
#include <stdio.h>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>  */

//using namespace std;
//#include "../icenReg_files/basicUtilities.cpp"


class node_info{
public:
    vector<int> l;      //vector that indicates the observations for which this node is the left side
    vector<int> r;      //vector that indicated the observations for which this node is the right side
//    double par;         //log cumulative hazard
};

class obInf{
public:
    int l,r;
    double pob;
};


class icm_Abst{
public:
    void update_p_ob(int i);    //done, not checked
    
    double sum_llk();           //done, not checked
    // calculates the entire likelihood function.
    // Does not update eta or hazards!
    
    double par_llk(int ind);     //done, not checked
    // only calculates partial likelihood based on an active index
    
    vector<obInf> obs_inf;
    vector<node_info> node_inf;
    
    void numericBaseDervsAllRaw(vector<double> &d1, vector<double> &d2);
    
    void icm_addPar(vector<double> &delta);

    void numericBaseDervsOne(int raw_ind, vector<double> &d);
    void numericBaseDervsAllAct(vector<double> &d1, vector<double> &d2);

    
    void update_etas();
	virtual void stablizeBCH() = 0;
    void recenterBCH();
	
    void icm_step();
    
    void numericRegDervs();
    void covar_nr_step();
    
    virtual double basHaz2CondS(double ch, double eta) = 0;     //done
    virtual double baseS2CondS(double s, double eta) = 0;
    virtual double base_d1_contr(double h, double pob, double eta) = 0; //done, not checked
    virtual double reg_d1_lnk(double ch, double xb, double log_p) = 0;
    virtual double reg_d2_lnk(double ch, double xb, double log_p) = 0;
    
    void calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1);
    void rawDervs2ActDervs();
    
    Eigen::VectorXd     baseCH;     //Vector of baseline log cumulative hazards.
                                    //baseH[0] fixed to -Inf, baseH[k-1] = Inf
	double intercept;				//used for numerical stabilization
	
    Eigen::VectorXd     backupCH;   //used to save values in optimization steps
    Eigen::VectorXd     propVec;    //used for proposition step during NR update on regression parameters
 /*   Eigen::VectorXd     H_d1;       //Vector of derivatives for CH's
    Eigen::MatrixXd     H_d2;       //Hessian for CH's          */
    Eigen::VectorXd     base_p_obs; //Baseline probability of each observation  //initialized
    Eigen::VectorXd     etas;       //linear combination of regression parameters   //initialized
    Eigen::VectorXd     expEtas;    //exp(etas) //initialized
    Eigen::VectorXd     reg_par;    //regression parameters //initialized
    Eigen::MatrixXd     covars;     //covariates        //initialized
    Eigen::VectorXd     reg_d1;     //first derivatives of regression parameters        //initialized
    Eigen::MatrixXd     reg_d2;     //Hessian for derivatives       //initialized
//    Eigen::VectorXd     reg_d2;     //second derivatives: ignoring off diagonals!
    
    vector<double> w;
    
    double maxBaseChg;      //Max change in baseline parameters during icm step
    double h;
    bool hasCovars;
    bool updateCovars;
    
    bool startGD;
    vector<double> baseS;
    vector<double> baseP;
    vector<double> baseP_backup;
    vector<double> d_cond_S_left;
    vector<double> d_cond_S_right;
    vector<double> base_p_derv;
    vector<double> base_p_derv2;			// For computing 2nd derivative
    vector<double> base_p_2ndDerv;
    vector<double> prop_p;
    double llk_from_p();
    double numeric_p_der(int i);
    
    double dervConS_fromBaseS(double s, double eta);
    void baseCH_2_baseS();
    void baseS_2_baseP();
    void baseP_2_baseS();
    void baseS_2_baseCH();
    void calc_cond_S_derv();
    void calc_base_p_derv();
    double getMaxScaleSize( vector<double> &p, vector<double> &prop_p);
    void gradientDescent_step();
    void experimental_step();
    void EM_step();
    
    vector<double> dob_dp_both;
    vector<double> dob_dp_rightOnly;

	double run(int maxIter, double tol, bool useGA, int baselineUpdates);
    
    void numeric_dobs_dp(bool forGA);
    void numeric_dobs2_d2p();
    
    double cal_log_obs(double s1, double s2, double eta);
    
    vector<bool> usedVec;
    
    double almost_inf;
    int failedGA_counts;
    int iter;
    int numBaselineIts;
    bool useFullHess;
    
    double exchangeAndUpdate(double delta, int i1, int i2);
    // REQUIRES baseP BEING UP TO DATE!!!
    
    vector<int> exchangeIndices;
    
    void checkCH();
    
    void last_p_update();
    void vem();
    void exchange_p_opt(int i1, int i2);
    void vem_sweep();
    void vem_sweep2();
};

void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, icm_Abst* icm_obj);
//function for setting up a actSet_Abst class

void cumhaz2p_hat(Eigen::VectorXd &ch, vector<double> &p);

class icm_ph : public icm_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        return(exp(-exp(ch + eta) )) ;}
    
    double baseS2CondS(double s, double eta){
        if(s >= 1.0) return(1.0);
        if(s <= 0.0) return(0.0);
/*        double expEta = exp(eta);
        double ans = pow(s, expEta);    */
        double logCH = log( -log(s) );
        double ans = exp(-exp(logCH + eta));
        return(ans);
    }
    
    double base_d1_contr(double ch, double pob, double eta){
        double expVal = -exp(eta + ch);
        double logAns = eta + ch + expVal - log(pob);
        return (-exp(logAns));
    }
    
    double reg_d1_lnk(double ch, double xb, double log_p){
        double term1 = -exp(ch + xb);
        return(-exp(term1 + ch + xb - log_p));
    }
    double reg_d2_lnk(double ch, double xb, double log_p){
        double term1 = -exp(ch + xb);
        double term2 = exp(term1 - log_p);
        return(term1 * term2 + term1 * term1 *term2);
    }
	
	void stablizeBCH(){
		int k = baseCH.size();
		double thisChange = baseCH[k-2] - 2.0;
		intercept += thisChange;
		for(int i = 1; i < (k-1); i++){
			baseCH[i] -= thisChange;
		}
		update_etas();	
	}
	
	
    virtual ~icm_ph(){};
};


class icm_po : public icm_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        double mu = exp(ch);
        double s = exp(-mu);
        double s_nu = exp(eta - mu);
        return( (s_nu) / (s_nu - s + 1)) ;}
    
    double baseS2CondS(double s, double eta){
        double nu = exp(eta);
        double s_nu = s * nu;
        return((s_nu)/ (s_nu - s + 1) );
    }
    
    double base_d1_contr(double ch, double pob, double eta){
        double s = exp(-exp(ch));
        double s_nu = exp(eta - exp(ch));
        double logAns = -log(pob) - 2 * log(s_nu - s + 1) + ch - exp(h);
        return (-exp(logAns));
    }
    
    double reg_d1_lnk(double ch, double xb, double log_p){
        double s = exp(-exp(ch));
        double a = exp(xb-exp(ch));
        double ans = exp( log(a *(1-s)) - 2 * log( a - s + 1) - log_p);
        return(ans);
//        return( a * (1-s) / pow(a - s + 1, 2.0) );
    }
    double reg_d2_lnk(double ch, double xb, double log_p){
        double s = exp(-exp(ch));
        double a = exp(xb-exp(ch));
        double b = a - s + 1;
        double top =  (a * (1 - s) * b - 2 * a * a *(1-s)) ;
        double bottom = pow(b, 3.0);
        double ans = top/(bottom * exp(log_p));
        return(ans);
    }
	void stablizeBCH(){}
	
    virtual ~icm_po(){};
};

extern "C" {
SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType,
 			  SEXP R_w, SEXP R_use_GD, SEXP R_maxiter,
 			  SEXP R_baselineUpdates, SEXP R_useFullHess, SEXP R_updateCovars,
 			  SEXP R_initialRegVals);
    SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals);
}
#endif /* defined(____ic_sp_cm__) */
