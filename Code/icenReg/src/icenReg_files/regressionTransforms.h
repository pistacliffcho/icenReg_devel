//
//  regressionTransforms.h
//  
//
//  Created by Cliff Anderson Bergman on 11/17/15.
//
//

#ifndef ____regressionTransforms__
#define ____regressionTransforms__


class condProbCal{
public:
    SEXP* baselineInfo;
    double (*getBaseSurv)(double q, SEXP bli);
    double (*transformSurv)(double s, double nu);
    double (*getBaseQ)(double s, SEXP bli);
    double (*transform_p)(double q, double nu);
    vector<double> preppedParams;
    condProbCal(SEXP regType, SEXP baseType, SEXP bli);
    bool isOK;
};


// SURVIVAL TRANSFORMS
double propHazTrans(double s, double nu){
    double ans = pow(s, nu);
    return( ans );
}

double propOddsTrans(double s, double nu){
    double ans;
    double prod = s * nu;
    ans = prod/(prod - s + 1);
    return(ans);
}

double transform_p_ph(double p, double nu){
    double log_s = log(1.0 - p);
    double log_trans = log_s / nu;
    return(1.0 - exp(log_trans));
}
double transform_p_po(double p, double nu){
    double s = 1.0 - p;
    return(1.0 - s * (1/nu) / (s * 1/nu - s +1));
}



// BASELINE MODELS
double getExpSurv(double q, SEXP bli){
    return( pexp(q, exp(REAL(bli)[0]), 0, 0) );
}
double getWeibSurv(double q, SEXP bli){
    return( pweibull(q, exp(REAL(bli)[0]), exp(REAL(bli)[1]), 0, 0));
}
double getLogNormSurv(double q, SEXP bli){
    return( pnorm(log(q), REAL(bli)[0], exp(REAL(bli)[1]), 0, 0) );
}
double getGammaSurv(double q, SEXP bli){
    return( pgamma(q, exp(REAL(bli)[0]), exp(REAL(bli)[1]), 0,0 ) );
}
double getLgLgsticSurv(double q, SEXP bli){
    return( 1.0 - ic_ploglogistic(q, exp(REAL(bli)[0]), exp(REAL(bli)[1])));
}
double getGenGammaSurv(double q, SEXP bli){
    return(1.0 - ic_pgeneralgamma(q, REAL(bli)[0], exp(REAL(bli)[1]), REAL(bli)[2]));
}
double getNonParSurv(double q, SEXP bli);




double getExpQ(double p, SEXP bli){
    return( qexp(p, exp(REAL(bli)[0]), 1, 0) );
}

double getWeibQ(double p, SEXP bli){
    return( qweibull(p, exp(REAL(bli)[0]), exp(REAL(bli)[1]), 1, 0));
}
double getLogNormQ(double p, SEXP bli){
    return( exp( qnorm(p, REAL(bli)[0], exp(REAL(bli)[1]), 1, 0) ) );
}
double getGammaQ(double p, SEXP bli){
    return( qgamma(p, exp(REAL(bli)[0]), exp(REAL(bli)[1]), 1,0 ) );
}
double getLgLgsticQ(double p, SEXP bli){
    return(ic_qloglogistic(p, exp(REAL(bli)[0]), exp(REAL(bli)[1]) ) );
}
double getGenGammaQ(double p, SEXP bli){
    return(ic_qgeneralgamma(p, REAL(bli)[0], exp(REAL(bli)[1]), REAL(bli)[2] ) );
}

double getNonParQ(double q, SEXP bli);



extern "C" {
    SEXP s_regTrans(SEXP times, SEXP etas,
                    SEXP bli, SEXP regType, SEXP baseType);
    SEXP q_regTrans(SEXP q, SEXP etas,
                    SEXP bli, SEXP regType, SEXP baseType);
}
#endif /* defined(____regressionTransforms__) */
