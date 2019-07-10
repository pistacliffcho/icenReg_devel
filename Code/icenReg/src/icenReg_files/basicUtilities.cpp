//
//  basicUtilities.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//

#include "myPAVAalgorithm.cpp"
#define SMALLNUMBER 0.00000000001;

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

void add_vec(double a, vector<double> &vec){
    int thisSize = vec.size();
    for(int i = 0; i < thisSize; i++)
        vec[i] += a;
}

void add_vec(vector<double> &a, vector<double> &vec){
    int thisSize = vec.size();
    int thisSize2 = a.size();
    if(thisSize != thisSize2){
        Rprintf("warning: sizes do not match in add_vec\n");
        return;
    }
    for(int i = 0; i < thisSize; i++)
        vec[i] += a[i];
}

void add_vec(double lambda, vector<double> &a, vector<double> &vec){
    int thisSize = vec.size();
    int thisSize2 = a.size();
    if(thisSize != thisSize2){
        Rprintf("warning: sizes do not match in add_vec\n");
        return;
    }
    for(int i = 0; i < thisSize; i++)
        vec[i] += a[i] * lambda;
}


double signVal(double x){
    if(x > 0) return 1.0;
    return -1.0;
}

void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat){
    SEXP Rdims = Rf_getAttrib(r_mat, R_DimSymbol);
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

Rcpp::NumericMatrix eigen2RMat(Eigen::MatrixXd &e_mat){
	int nRow = e_mat.rows();
	int nCol = e_mat.cols();
	Rcpp::NumericMatrix ans(nRow, nCol);
	for(int i = 0; i < nRow; i++){
		for(int j = 0; j < nCol; j++){
			ans(i,j) = e_mat(i,j);
		}
	}
	return(ans);
}

void Rvec2eigen(SEXP r_vec, Eigen::VectorXd &e_vec){
    int k = LENGTH(r_vec);
    e_vec.resize(k);
    for(int i = 0; i < k; i++)
        e_vec[i] = REAL(r_vec)[i];
}

Rcpp::NumericVector eigen2RVec(Eigen::VectorXd &e_vec){
	int k = e_vec.size();
	Rcpp::NumericVector ans(k);
	double* cptr = &ans[0];
	for(int i = 0; i < k; i++){
		cptr[i] = e_vec[i];
	}
	return(ans);
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

double ic_qloglogistic(double p, double a, double b){
    double expVal = -1.0/b;
    double odds = (1.0 - p)/p;
    double ans = a * pow(odds, expVal);
    return(ans);
}

double ic_dlnorm(double x, double mu, double s){
    double denom = x * s * sqrt(2 * M_PI);
    double expPart = pow((log(x) - mu), 2) / (2 * s * s);
    double ans = exp(-expPart) / denom;
    return (ans);
}

double ic_plnorm(double x, double mu, double s){
    return(R::pnorm(log(x), mu, s, 0, 0));
}




double ic_dgeneralgamma(double x, double mu, double s, double Q){
    if(Q == 0){
        return(ic_dlnorm(x, mu, s));
    }
    double y = log(x);
    double w = ((y - mu)/s);
    double Q_inv = 1/(Q * Q);
    double logAns = -log(s*x) + log(abs(Q)) + Q_inv * log(Q_inv) + Q_inv * (Q*w - exp(Q*w)) - Rf_lgammafn(Q_inv);
    return(exp(logAns));
}

double ic_pgeneralgamma(double q, double mu, double s, double Q){
    if(Q == 0){ return(ic_plnorm(q, mu, s));}
    double y = log(q);
    double w = ((y-mu)/s);
    double Q_inv = 1 / (Q*Q);
    double expnu = exp(Q*w) * Q_inv;
    double ans;
    if(Q > 0){ ans = 1 - Rf_pgamma(expnu, Q_inv, 1.0,  0, 0);}
    else{ ans = Rf_pgamma(expnu, Q_inv, 1.0, 0, 0);}
    return(ans);
}

double ic_qgeneralgamma(double p, double mu, double s, double Q){
    if(Q == 0){ return(R::qlnorm(p, mu, s, 1, 0)); }
    double Q2 = (Q*Q);
    double Q_inv = 1/Q2;
    double part2 = s * (log(Q2 * R::qgamma(p, Q_inv, 1.0, 0, 0)) / Q);
    double ans = exp(mu + part2);
    return(ans);
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


void pavaForOptim(Eigen::VectorXd &d1, Eigen::VectorXd  &d2, Eigen::VectorXd  &x, Eigen::VectorXd  &prop_delta){
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
    if(r == max){ 
    	points.push_back(r);
    	return;
    }
    int thisSize = points.size();
    for(int i = 0; i < thisSize; i++){
        if(l < points[i] && (r + 1) >= points[i]) { 
        	return;
        }
    }
	points.push_back(r+1);
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
        
        SEXP ans = Rf_allocVector(REALSXP, k);
        PROTECT(ans);
        for(int i = 0; i<k; i++){ REAL(ans)[i] = prop_delta[i];}
        UNPROTECT(1);
        return(ans);
    }
}



extern "C"{
    SEXP fastNumericInsert(SEXP newVals, SEXP target, SEXP indices){
        int k1 = LENGTH(newVals);
        int k2 = LENGTH(indices);
        int t_k = LENGTH(target);
        
        int* c_inds = INTEGER(indices);
        double* c_targ = REAL(target);
        double* nv = REAL(newVals);
        
        for(int i = 0; i < k2; i++){
            if(c_inds[i] > t_k){
                Rprintf("error: index too large for target in fastNumericInsert\n");
                return(target);
            }
            if(c_inds[i] < 1){
                Rprintf("error: index less than 1 for target in fastNumericInsert\n");
                return(target);
            }
        }
        
        if(k1 == 1){
            double onlyVal = nv[0];
            for(int i = 0; i < k2; i++){
                c_targ[c_inds[i] - 1] = onlyVal;
            }
        }
        else{
            for(int i = 0; i < k2; i++){
                c_targ[c_inds[i] - 1] = nv[i];
                }
            }
    return(target
           );
    }
}


//Added for bivariate NPMLE

void SEXP2doubleVec(SEXP R_vec, std::vector<double> &c_vec){
    int k = LENGTH(R_vec);
    c_vec.resize(k);
    for(int i = 0; i < k; i++) { c_vec[i] = REAL(R_vec)[i]; }
}

void doubleVec2SEXP(std::vector<double> &c_vec, SEXP R_vec){
    int k = c_vec.size();
    int r_k = LENGTH(R_vec);
    if(k != r_k){ Rprintf("Warning: doubleVec2SEXP skipped due to different lengths\n");}
    else {for(int i = 0; i < k; i++){ REAL(R_vec)[i] = c_vec[i];} }
}

void SEXPIndex2intIndex(SEXP R_Inds, std::vector<int> &c_inds){
    int k = LENGTH(R_Inds);
    c_inds.resize(k);
    for(int i = 0; i < k; i++){
        c_inds[i] = INTEGER(R_Inds)[i] - 1;
    }
}

void indexVec2SEXP(std::vector<int> &c_vec, SEXP R_vec){
    int k = c_vec.size();
    int r_k = LENGTH(R_vec);
    if(k != r_k){ Rprintf("Warning: doubleVec2SEXP skipped due to different lengths\n");}
    else {for(int i = 0; i < k; i++){ INTEGER(R_vec)[i] = c_vec[i] + 1;} }
}


void findIndexDiffs(std::vector<int> &in1, std::vector<int> &in2,
                    std::vector<int> &in1not2, std::vector<int> &in2not1){
    int k1 = in1.size();
    int k2 = in2.size();
    
    in1not2.resize(0);
    in2not1.resize(0);
    
    if(k1 == 0 || k2 == 0) {return;}
    
    in1not2.reserve(k1);
    in2not1.reserve(k2);
    
    int i1 = 0;
    int i2 = 0;
    int ind1;
    int ind2 = in2[i2];
    
    while(i1 < k1){
        ind1 = in1[i1];
        while(ind2 < ind1){
            in2not1.push_back(ind2);
            i2++;
            if(i2 < k2){ ind2 = in2[i2]; }
            else {ind2 = in1[k1 -1] + 1;}
        }
        if(ind1 < ind2) { in1not2.push_back(ind1); }
        else{
            i2++;
            if(i2 < k2){ ind2 = in2[i2]; }
            else {ind2 = in1[in1.size() -1] + 1;}
        }
        i1++;
    }
    
    while(i2 < k2){
        ind2 = in2[i2];
        in2not1.push_back(ind2);
        i2++;
    }
    /*
     int difk1 = in1not2.size();
     int difk2 = in2not1.size();
     if(difk1 == 0 || difk2 == 0){
     Rprintf("warning: difk1 = %d, difk2 = %d\n", difk1, difk2);
     for(int i = 0; i < k1; i++){Rprintf("%d ", in1[i]);}
     Rprintf("\n");
     for(int i = 0; i < k2; i++){Rprintf("%d ", in2[i]);}
     Rprintf("\n");
     }   */
}



void drop_index(int d_ind, std::vector<int> &indVec){
    int k = indVec.size();
    for(int i = 0; i < k; i++){
        if(d_ind == indVec[i]){
            indVec.erase( indVec.begin() + i);
            return;
        }
    }
    Rprintf("error: attempting to drop an index not found in vector\n");
}

void add_index(int a_ind, std::vector<int> &indVec){
    int k = indVec.size();
    if(k == 0){
        indVec.insert(indVec.begin(), a_ind);
        return;
    }
    if(a_ind < indVec[0]){
        indVec.insert(indVec.begin(), a_ind);
        return;
    }
    for(int i = 1; i < k; i++){
        if(a_ind < indVec[i]){
            indVec.insert( indVec.begin() + i, a_ind);
            return;
        }
    }
    if(a_ind > indVec[indVec.size() -1] ){
        indVec.push_back(a_ind);
        return;
    }
    Rprintf("error: trying to insert index that is already in vector. Index = %d\n", a_ind);
}

std::vector<int> getSEXP_MatDims(SEXP R_mat){
    SEXP Rdims = Rf_getAttrib(R_mat, R_DimSymbol);
    PROTECT(Rdims);
    std::vector<int> ans(2);
    ans[0] = INTEGER(Rdims)[0];
    ans[1] = INTEGER(Rdims)[1];
    UNPROTECT(1);
    return(ans);
}

void getPosNegIndices(std::vector<double> &vals, std::vector<int> &isPos, std::vector<int> &isNeg){
    isPos.resize(0);
    isNeg.resize(0);
    int k = vals.size();
    for(int i = 0; i < k; i++){
        if(vals[i] > 0.0){isPos.push_back(i);}
        else{isNeg.push_back(i);}
    }
}

void getRelValIndices(double relVal, std::vector<double> &vals, std::vector<int> &subIndex,
                      std::vector<int> &above, std::vector<int> &below,
                      int *max, int *min){
    above.resize(0);
    below.resize(0);
    int k = vals.size();
    int k2 = subIndex.size();
    double maxVal = R_NegInf;
    double minVal = R_PosInf;
    
    (*max) = -1;
    (*min) = -1;
    
    if(k != k2){
        Rprintf("in getPosNegIndices, k != k2! Quiting.\n");
        return;
    }
    
    for(int i = 0; i < k; i++){
        if(vals[i] > relVal){
            above.push_back(subIndex[i]);
            if(maxVal < vals[i]){
                maxVal = vals[i];
                (*max) = subIndex[i];
            }
        }
        else{
            below.push_back(subIndex[i]);
            if(minVal > vals[i]){
                minVal = vals[i];
                (*min) = subIndex[i];
            }
        }
    }
    
}

double directional_derv(vector<double> &derv, vector<double> &delta){
    int k = derv.size();
    int k2 = delta.size();
    if(k != k2){
        Rprintf("warning: sizes don't match in directional_derv\n");
        return(0.0);
    }
     double ans = 0.0;
    for(int i = 0; i < k; i++){
        ans += derv[i] * delta[i];
    }
    return(ans);
}


void makeUnitVector(vector<double> &v){
    double sum = 0;
    int k = v.size();
    for(int i = 0; i < k; i++){
        sum += abs(v[i]);
    }
    for(int i = 0; i < k; i++){
        v[i] = v[i]/sum;
    }
}

void getUniqInts(int i1, int i2, vector<int> &uniqInts, vector<vector<int> > &vec_vec, vector<bool> &usedVec){
    uniqInts.clear();
    int tot_pos_lng = 0;
    int lv = vec_vec.size();
    if(i1 >= lv){Rprintf("i1 too big in getUniqInts\n");return;}
    if(i2 >= lv){Rprintf("i2 too big in getUniqInts\n");return;}
    
    for(int i = i1; i <= i2; i++){
        tot_pos_lng += vec_vec[i].size();
    }
    
    uniqInts.reserve(tot_pos_lng);
    vector<int>* theseIndices;
    int thisSize;
    int thisIndex;
    for(int i = i1; i <= i2; i++){
        theseIndices = &vec_vec[i];
        thisSize = theseIndices->size();
        for(int j = 0; j < thisSize; j++){
            thisIndex = (*theseIndices)[j];
            if(usedVec[thisIndex] == false){
                usedVec[thisIndex] = true;
                uniqInts.push_back(thisIndex);
            }
        }
    }
    thisSize = uniqInts.size();
    for(int i = 0; i < thisSize; i++){
        usedVec[i] = false;
    }
}

int isValueInInterval(double val, double l, double r){
	if(val < l) return(-1);
	if(val > r) return(1);
	return(0);
}

int isValueInInterval(double val, int ind, 
					  vector<double>& lvec, vector<double>& rvec){
	return(isValueInInterval(val, lvec[ind], rvec[ind]));					  
}

int findSurroundingVals(double val, vector<double>& leftVec,
						vector<double>& rightVec, bool isLeft){
	
	int a = 0;
	int b = leftVec.size()-1;
	if(b == 0){return(0);}
	if(isValueInInterval(val, R_NegInf, rightVec[0]) == 0) return(0);
	if(isValueInInterval(val, leftVec[b], R_PosInf) == 0) return(b);
	
/*	a++;
	b--;	*/
	
	int maxTries = b;
	
	int propInd = (a + b)/2;
	int tries = 0;
	int testVal;
	while( b - a > 1 && tries < maxTries){
		tries++;
		propInd = (a + b)/2;
		testVal = isValueInInterval(val, propInd, leftVec, rightVec);
		if(testVal == 0){ return(propInd);}
		if(testVal == -1){ b = propInd; }
		else{ a = propInd; }
	}
	if(a == b){
		Rprintf("this is very surprising... a = %d, size = %d\n", a, leftVec.size());
		return(a);
	}
	if( isLeft ) return(b); 
	return(a);
}

Eigen::MatrixXd xtx(Eigen::MatrixXd &x){
	int nCols = x.cols();
	int nRows = x.rows();
	Eigen::MatrixXd ans(nCols, nCols);
	for(int i = 0; i < nCols; i++){
		for(int j = 0; j <= i; j++){
			ans(i,j) = 0;
			for(int r = 0; r < nRows; r++){ ans(i,j) += x(r, i) * x(r, j);}
			ans(j,i) = ans(i,j);
		}
	}
	return(ans);
}

Eigen::MatrixXd xtx(Eigen::MatrixXd &x, int row_start, int row_end){
	int nCols = x.cols();
	Eigen::MatrixXd ans(nCols, nCols);
	for(int i = 0; i < nCols; i++){
		for(int j = 0; j <= i; j++){
			ans(i,j) = 0;
			for(int r = row_start; r < row_end; r++){ ans(i,j) += x(r, i) * x(r, j);}
			ans(j,i) = ans(i,j);
		}
	}
	return(ans);
}

Eigen::MatrixXd copyRows(Eigen::MatrixXd &x, int row_start, int row_end){
	int nCols = x.cols();
	int nRows = row_end - row_start + 1;
	
	Eigen::MatrixXd ans(nRows, nCols);
	int cur_row = 0;
	for(int i = row_start; i <= row_end; i++){
		ans.row(cur_row) = x.row(i);
		cur_row++;
	}
	return(ans);
}
