//
//  bivariateNPMLE.h
//  
//
//  Created by Cliff Anderson Bergman on 8/29/15.
//
//

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <vector>

extern "C"{
    SEXP R_testDiffStep(SEXP in1, SEXP in2);
    SEXP optCliq(SEXP cliqMat, SEXP R_tol, SEXP R_innerLoops, SEXP R_outerLoops, SEXP R_obs_as_rows);
}

/*
 void squeezeStep(double *p1, double *p2,
                 std::vector<int> &in1Not2, std::vector<int> &in2Not1,
                 std::vector<double> &p_obs);
*/

void SEXPMat2pmass_info(SEXP r_mat, std::vector< std::vector<int> > &cliqMat, bool obs_as_rows );
double compute_pob(std::vector<int> &this_in, std::vector<double> &p_vec);

//void SEXPMat2cliqMat(SEXP r_mat, std::vector< std::vector<int> > cliqMat, bool obs_as_rows);

class bvcen {
public:
    std::vector< std::vector<int> > pmass_in_ob;
//    std::vector< std::vector<int> > org_pmass_in_ob;
    //each vector should contain the information about a given pmass
    
    
/*    std::vector< std::vector<int> > pob_in_pmass;
    //each vector should contain the information about a given pobs
    
    std::vector< std::vector<int> > org_pob_in_pmass;
    //full information      */
    
    std::vector<int> pos_pmass;
    std::vector<int> zer_pmass;
    std::vector<double> p_mass;
    std::vector<double> p_obs;
    
    std::vector<int> drop_inds;
    std::vector<double> p_obs_inv;
    
    std::vector<double> dp_full;
    std::vector<double> dp_act;
    void calc_full_dp();
    void calc_act_dp();
    void full_em();
    void act_em();

    std::vector<int> negInds;
    std::vector<int> posInds;
    
    double fullError;
    double actError;
    
    void drop_pmass(int p_ind);
    void add_pmass(int p_ind);

    void drop_zeros();
    void add_points();
    
    
    double llk();
    void update_pobs();
    void squeeze(int i, int j);
    void squeezeInternal(int p1_ind, int p2_ind, std::vector<int> &in1not2, std::vector<int> &in2not1);
    
    void vem_act();
    
    double get_ptot();
    
    bvcen(SEXP cliqMat, SEXP R_obs_as_rows);
};
