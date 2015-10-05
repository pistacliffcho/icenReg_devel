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

/*class actPointInf{
public:
    int ind;                   //beginning and ending of parameters affected by adjusting active points
    double par;
    vector<int> dep_obs;
//    vector<int> dep_nodes;      //don't think this is necessary; l, r tells us all this!
};  */

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
    
    void icm_step();
    
    
    void numericRegDervs();
    void covar_nr_step();
    
     
    virtual double basHaz2CondS(double ch, double eta) = 0;     //done
    virtual double base_d1_contr(double h, double pob, double eta) = 0; //done, not checked
    virtual double reg_d1_lnk(double ch, double xb, double log_p) = 0;
    virtual double reg_d2_lnk(double ch, double xb, double log_p) = 0;
    
    void calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1);
    void rawDervs2ActDervs();
    
    Eigen::VectorXd     baseCH;     //Vector of baseline log cumulative hazards.
                                    //baseH[0] fixed to -Inf, baseH[k-1] = Inf
    Eigen::VectorXd     H_d1;       //Vector of derivatives for CH's
    Eigen::MatrixXd     H_d2;       //Hessian for CH's
    Eigen::VectorXd     base_p_obs; //Baseline probability of each observation  //initialized
    Eigen::VectorXd     etas;       //linear combination of regression parameters   //initialized
    Eigen::VectorXd     expEtas;    //exp(etas) //initialized
    Eigen::VectorXd     reg_par;    //regression parameters //initialized
    Eigen::MatrixXd     covars;     //covariates        //initialized
    Eigen::VectorXd     reg_d1;     //first derivatives of regression parameters        //initialized
    Eigen::MatrixXd     reg_d2;     //Hessian for derivatives       //initialized

    vector<double> w;
    
    double maxBaseChg;      //Max change in baseline parameters during icm step
    double h;
    bool hasCovars;
};

void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, icm_Abst* icm_obj);
//function for setting up a actSet_Abst class

//void addDepNodes(vector<int> &intoVec, int l, int r, vector<node_info> &nf);
//

void cumhaz2p_hat(Eigen::VectorXd &ch, vector<double> &p);

class icm_ph : public icm_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        return(exp(-exp(ch + eta) )) ;}
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
    virtual ~icm_po(){};
};

extern "C" {
    SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType, SEXP R_w);
    SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals);
}
#endif /* defined(____ic_sp_cm__) */
