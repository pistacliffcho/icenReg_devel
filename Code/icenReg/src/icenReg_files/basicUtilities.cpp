//
//  basicUtilities.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#include "myPAVAalgorithm.cpp"

double max(double a, double b){
    if(a > b) return(a);
    return(b);
}

double max(int a, int b){
    if(a > b) return(a);
    return(b);
}

template <typename T>
int getMaxIndex(T &v){
    double mVal = R_NegInf;
    int mInd = 0;
    int k = v.size();
    for(int i = 0; i < k; i++){
        if(v[i] == max(mVal, v[i])) {
            mVal = v[i];
            mInd = i;
        }
    }
    return(mInd);
}

double min(double a, double b){
    if(a < b) return (a);
    return(b);
}

double max(Eigen::VectorXd v){
    double maxVal = R_NegInf;
    for(int i = 0; i < v.size(); i++) {maxVal = max(maxVal, v[i]);}
    return(maxVal);
}

void mult_vec(double a, vector<double> &vec){
    int thisSize = vec.size();
    for(int i = 0; i < thisSize; i++)
        vec[i] *= a;
}

double signVal(double x){
    if(x > 0) return 1.0;
    return -1.0;
}

void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat){
    SEXP Rdims = getAttrib(r_mat, R_DimSymbol);
    PROTECT(Rdims);
    int nRows = INTEGER(Rdims)[0];
    int nCols = INTEGER(Rdims)[1];
    
    e_mat.resize(nRows, nCols);
    for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nCols; j++)
            e_mat(i, j) = REAL(r_mat)[i + j*nRows];
    }
    UNPROTECT(1);
}

void Rvec2eigen(SEXP r_vec, Eigen::VectorXd &e_vec){
    int k = LENGTH(r_vec);
    e_vec.resize(k);
    for(int i = 0; i < k; i++)
        e_vec[i] = REAL(r_vec)[i];
}



double ic_dloglogistic(double x, double a, double b){
    double x_a = x/a;
    double x_a_b = pow(x_a, b);
    double ans = (b/a) * (x_a_b/x_a) / pow(1 + x_a_b, 2);
    return(ans);
}

double ic_ploglogistic(double x, double a, double b){
    return( 1/ (1 + pow(x/a, -b)));
}

double ic_dlnorm(double x, double mu, double s){
    double denom = x * s * sqrt(2 * PI);
    double expPart = pow((log(x) - mu), 2) / (2 * s * s);
    double ans = exp(-expPart) / denom;
    return (ans);
}

double ic_plnorm(double x, double mu, double s){
    return(pnorm(log(x), mu, s, 0, 0));
}


void pavaForOptim(vector<double> &d1, vector<double> &d2, vector<double> &x, vector<double> &prop_delta){
    int k = d1.size();
    int d2_size = d2.size();
    int x_size = x.size();
    if(k != d2_size || k!= x_size){ Rprintf("incorrect sizes provided to pavaForOptim\n"); return;}
    prop_delta.resize(k);
    vector<double> y(k);
    vector<double> w(k);
    
    for(int i = 0; i < k; i++){
        y[i] = -d1[i]/d2[i] + x[i];
        w[i] = d2[i]/2;
    }
    int k_sign = k;
    pava( &y[0], &w[0], &k_sign );
    for(int i = 0; i < k; i++){
        prop_delta[i] = y[i] - x[i];
    }
}

void addIfNeeded(vector<int> &points, int l, int r, int max){
    if(r > max)     {Rprintf("warning: r > max\n"); return;}
    bool chg_btw = false;
    int thisSize = points.size();
    for(int i = 0; i < thisSize; i++){
        if(l < points[i] && (r + 1) >= points[i]) { chg_btw = true;};
 //       if(l <= points[i])  l_below = true;
 //       if( r >= points[i]) r_geq   = true;
    }
    if(chg_btw) {return;}   //no need to add point
    if(r == max)    points.push_back(r);
    else            points.push_back(r+1);
}




extern "C"{
    SEXP pava(SEXP R_d1, SEXP R_d2, SEXP R_x){
        int k = LENGTH(R_d1);
        if(k!= LENGTH(R_d2) || k != LENGTH(R_x) ){Rprintf("sizes don't match! Quiting pava\n"); return(R_NilValue);}
        vector<double> d1, d2, x, prop_delta;
        d1.resize(k);
        d2.resize(k);
        x.resize(k);
        for(int i = 0; i<k; i++){
            d1[i] = REAL(R_d1)[i];
            d2[i] = REAL(R_d2)[i];
            x[i]  = REAL(R_x)[i];
        }
        pavaForOptim(d1, d2, x, prop_delta);
        
        SEXP ans = allocVector(REALSXP, k);
        PROTECT(ans);
        for(int i = 0; i<k; i++){ REAL(ans)[i] = prop_delta[i];}
        UNPROTECT(1);
        return(ans);
    }
}
