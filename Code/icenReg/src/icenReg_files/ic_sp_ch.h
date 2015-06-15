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

class actPointInf{
public:
    int ind;                   //beginning and ending of parameters affected by adjusting active points
    double par;
    vector<int> dep_obs;
//    vector<int> dep_nodes;      //don't think this is necessary; l, r tells us all this!
};

class actSet_Abst{
public:
    void update_p_ob(int i);    //done, not checked
    
    double sum_llk();           //done, not checked
    // calculates the entire likelihood function.
    // Does not update eta or hazards!
    
    double par_llk(int a_ind);     //done, not checked
    // only calculates partial likelihood based on an active index
    
    void act_setPar(int act_i, double val);     //done, not checked
    //updates the active index act_i according to an active set update
    //updates baseline hazards and dependent observations
    
    void act_addPar(int act_i, double delta){
        double newVal = actIndex[act_i].par + delta;
        act_setPar(act_i, newVal);
    }
    
    void act_addPar(vector<double> &delta);
    
    void addDepNodes(vector<int> &intoVec, int l, int r);
    // adds dependent obs into vector for use in par_llk
    
    vector<obInf> obs_inf;
    vector<node_info> node_inf;
    
    vector<actPointInf> actIndex;
    //vector of active point infos
    int getNextRawActInd(int act_i);
    
    void update_etas();

    
    void uniActiveOptim(int raw_ind);
    //univariate update of active points

    void numericBaseDervsOne(int raw_ind, vector<double> &d);
    void numericBaseDervsAllAct(vector<double> &d1, vector<double> &d2);
    void numericBaseDervsAllRaw(vector<double> &d1, vector<double> &d2);

    void analyticBase1stDerv(vector<double> &d1);
    
    void vem_step();
    void icm_step();
    
    
    void numericRegDervs();
    void covar_nr_step();
    
    int getActInd(int raw_ind);
    //takes in a raw index and gives back the corresponding active index.
    //returns -1 if not an index
    
    int getRawInd(int act_ind){return(actIndex[act_ind].ind);}
    //returns the raw index from the i_th active index
    
    
    
    void addActive(int raw_ind, double par);        //done
    //add active point

    void addActive(int raw_ind){ addActive(raw_ind, baseCH[raw_ind]);}
    // adds active point without adjusting it
    
    void removeActive(int act_ind); //done
    //remove active point
    
    int getNextActRawInd(int raw_ind);
    
    void checkIfActShouldDrop(int act_ind); //checks an active point and sees if it should be dropped
    
    virtual double basHaz2CondS(double ch, double eta) = 0;     //done
    virtual double base_d1_contr(double h, double pob, double eta) = 0; //done, not checked
    virtual double reg_d1_lnk(double ch, double xb) = 0;
    virtual double reg_d2_lnk(double ch, double xb) = 0;
    
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
    
    double h;
    bool hasCovars;
};

void setupActSet(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, actSet_Abst* actSet);
//function for setting up a actSet_Abst class

void addDepNodes(vector<int> &intoVec, int l, int r, vector<node_info> &nf);
//

void cumhaz2p_hat(Eigen::VectorXd &ch, vector<double> &p);

class actSet_ph : public actSet_Abst{
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
    
    double reg_d1_lnk(double ch, double xb){
        double term1 = -exp(ch + xb);
        return(-exp(term1 + ch + xb));
    }
    double reg_d2_lnk(double ch, double xb){
        double term1 = -exp(ch + xb);
        double term2 = exp(term1);
        return(term1 * term2 + term1 * term1 *term2);
    }
    virtual ~actSet_ph(){};
};


class actSet_po : public actSet_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        double s = exp(-exp(ch));
        double nu = exp(eta);
        return( (s*nu) / (nu * s - s + 1)) ;}
    double base_d1_contr(double ch, double pob, double eta){
        double s = exp(-exp(ch));
        double nu = exp(eta);
        double logAns = -log(pob) - 2 * log(s*nu - s + 1) + ch - exp(h);
        return (-exp(logAns));
    }
    
    double reg_d1_lnk(double ch, double xb){
        double s = exp(-exp(ch));
        double a = s * exp(xb);
        return( a * (1-s) / pow(a - s + 1, 2.0) );
    }
    double reg_d2_lnk(double ch, double xb){
        double s = exp(-exp(ch));
        double a = s * exp(xb);
        double b = a - s + 1;
        double top = a * (1 - s) * b - 2 * a * a *(1-s);
        double bottom = b * b * b;
        double ans = top / bottom;
        return(ans);
    }
    virtual ~actSet_po(){};
};

extern "C" {
    SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType, SEXP R_w);
    SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals);
}
#endif /* defined(____ic_sp_cm__) */
