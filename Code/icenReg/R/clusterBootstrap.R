# Generates raw indices for a single bootstrap resample
make_sample_inds = function(ID){
  unq_IDs = unique(ID)
  n_unq = length(unq_IDs)
  if(n_unq < 2){
    stop("IDs not unique")
  }
  raw_inds = seq_len(length(ID))
  ind_list = split(raw_inds, f = ID)
  ans = NULL
  id_sample = sample(1:n_unq, n_unq, replace = T)
  for(i in seq_len(n_unq) ){
    this_id = id_sample[i]
    ans = c(ans, ind_list[[this_id]])
  }
  return(ans)
}

# Compute the cluster bootstrap covariance
IR_cluster_vcov = function(fit, cluster_id, bs_samples = 1000){
  base_coefs = coef(fit)
  samps = matrix(nrow = bs_samples, 
                 ncol = length(base_coefs))
  
  fit_call = fit$call
  if(is(fit, "sp_fit"))
    fit_call$bs_samples = 0
  fit_call$data = as.name("BOOT__DATA")
  org_data = fit$getRawData()
  for(i in seq_len(bs_samples)){
    these_inds = make_sample_inds(cluster_id)
    BOOT__DATA = org_data[these_inds, ]
    this_fit = eval(fit_call)
    samps[i,] = coef(this_fit)
  }
  ans = cov(samps)
  colnames(ans) = names(base_coefs)
  rownames(ans) = names(base_coefs)
  return(ans)
}

#' Updates the covariance using cluster bootstrap
#' @param fit Either an ic_par or ic_sp model
#' @param ID Subject identifier
#' @param bs_samples Number of bootstrap samples
#' 
#' @description Adjusts error estimates for repeated measures data by use of the cluster bootstrap.
#' 
#' @details Standard models in icenReg assume independence between each observation. 
#' This assumption is broken if we can have multiple observations from a single subject, 
#' which can lead to an underestimation of the standard errors. \code{ir_clustBoot}
#' addresses this by using a cluster bootstrap to fix up the standard errors. 
#' 
#' Note that this requires refitting the model \code{bs_samples}, which means this can be 
#' fairly time consuming. 
#' 
#' @references 
#' Sherman, Michael, and Saskia le Cessie. 
#' "A comparison between bootstrap methods and 
#' generalized estimating equations for correlated 
#' outcomes in generalized linear models." 
#' Communications in Statistics-Simulation and Computation 26.3 (1997): 901-925.
#' 
#' @examples
#' # Simulating repeated measures data 
#' simdata = simIC_cluster(nIDs = 10, nPerID = 4)
#' 
#' # Fitting with basic model
#' fit = ic_par(cbind(l,u) ~ x1 + x2, data = simdata)
#' fit
#' 
#' # Updating covariance
#' ir_clustBoot(fit, ID = simdata$ID, bs_samples = 10)
#' # (Low number of bootstrap samples used for quick testing by CRAN, 
#' # never use this few!!)
#' 
#' # Note that the SE's have changed from above
#' fit
#' @export
ir_clustBoot = function(fit, ID, bs_samples = 1000){
  
  if(!(is(fit, "sp_fit") | is(fit, "par_fit") ) | is(fit, "ic_np")) 
     stop("fit must be output from ic_par or ic_sp")
  
  robust_vcov = IR_cluster_vcov(fit, ID, bs_samples)
  fit$var <- robust_vcov
  fit$llk = "NA: likelihood does not account for clustering"
  fit$depType = "Repeated Measures"
}

#' Simulates data with multiple observations per subject
#' @param nIDs Number of subjects
#' @param nPerID Number of observations per subject
#' @description Simulates data in which each subject is observed several times. 
#' In this case, the covariance matrix should be updated with \code{ir_clustBoot}.
#' @export
simIC_cluster = function(nIDs = 50, nPerID = 5){
  b1s = rnorm(nIDs, mean = 0.25, sd = 1)
  ans = NULL
  for(i in seq_len(nIDs)){
    this_data = simIC_weib(nPerID, b1 = b1s[i])
    this_data$ID = i
    ans = rbind(ans, this_data)
  }
  return(ans)
}