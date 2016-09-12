class aft_linkFun : public linkFun{
	// Note: this is a "fake" class!!
	// Because aft model behaves inherently differently than
	// prop odds and prop hazards models (does not alter survival prob,
	// but rather raw time), it does **not** need the same generic tools
	// as the other models. However, we will later build a class that does 
	// all the tools for this and still "looks" like an IC_parOpt object
	// So it will need a linkFun member...which only complains if it's 
	// ever called! 
	public:
	double con_s(double b_s, double nu){
		 Rprintf("Warning! con_s called from aft model.\n");
		 return(0.0);
	}
	double con_d(double b_d, double b_s, double nu){
		Rprintf("Warning! con_d called from aft model.\n");
		return(0.0);
	}
	double con_d_der(double b_d, double b_s, double nu){
		Rprintf("Warning! con_d_der called from aft model.\n");
		return(0.0);
	}
	double con_s_der(double b_s, double nu){
		Rprintf("Warning! con_s_der called from aft model.\n");
		return(0.0);
	}
	double cons_s_der2(double b_s, double nu){
		Rprintf("Warning! con_s_der2 called from aft model.\n");
		return(0.0);
	}
	virtual ~aft_linkFun(){};
};

class IC_parOpt_aft : public IC_parOpt{
public:
	double con_d(double base_t, double nu){
		double con_t = base_t / nu;
		double ans = blInf->base_d(con_t, b_pars);
		ans = ans / nu; 
		return(ans);
	}
	double con_s(double base_t, double nu){
		double con_t = base_t / nu;
		double ans = blInf->base_s(con_t, b_pars);
		return(ans);
	}
	
	void calculate_baseline_probs(){} 
	// baseline probs does not really fit into the aft
	// framework, so this is just a placeholder
		
	double calcLike_baseReady();
	void update_dobs_detas();
	virtual ~IC_parOpt_aft(){}
	IC_parOpt_aft(SEXP R_s_t, SEXP R_d_t, SEXP R_covars,
                     SEXP R_uncenInd, SEXP R_gicInd, SEXP R_lInd, SEXP R_rInd,
                     SEXP R_parType, SEXP R_linkType, SEXP R_w);
    IC_parOpt_aft(Rcpp::List R_list);
};	

