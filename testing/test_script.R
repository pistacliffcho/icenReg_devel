SAVE = FALSE
CHECK = TRUE
library(icenReg)
data(essIncData_small)
setwd('/Users/cliff/Desktop/icenReg_devel/testing')

makeTestObjects <- function(){
  set.seed(1)
  ans <- numeric()
  ans[1] <- ic_sp(cbind(inc_l, inc_u) ~ eduLevel * cntry, data = essIncData_small)$final_llk
  data('essIncData_small')
  ans[2] <- ic_sp(cbind(inc_l, inc_u) ~ eduLevel * cntry, data = essIncData_small, model = 'po')$final_llk
  ans[3] <- ic_par(cbind(inc_l, inc_u) ~ eduLevel * cntry, data = essIncData_small )$final_llk
  ans[4] <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, data = essIncData_small, model = 'po', dist = 'exponential')$final_llk
  ans[5] <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, data = essIncData_small, model = 'po', dist = 'weibull')$final_llk
  ans[6] <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, data = essIncData_small, model = 'po', dist = 'loglogistic')$final_llk
  ans[7] <- ic_par(Surv(inc_l, inc_u, type = 'interval2') ~ eduLevel * cntry, data = essIncData_small, model = 'po', dist = 'lnorm')$final_llk
  simdata <- simIC_weib(n = 200)
  fit <- ic_sp(cbind(l, u) ~ x1 + x2, data = simdata)
  ans[8] <- getFitEsts(fit)
  ans[9:10] <- getFitEsts(fit, p = c(.25, .75)) 
  newdata <- data.frame(x1 = c(1,2), x2 = c(-1,-2)) 
  ans[11:12] <- getFitEsts(fit, newdata, q = 1)
  ans_names <- c(rep('ic_sp_llk', 2), rep('ic_par_llk', 5), rep('fitted values', 5) )
  names(ans) <- ans_names
  return(ans)
}

if(SAVE){
  answers <- makeTestObjects()
  save(answers, file = 'answers.Rdata')
}

if(CHECK){
  load('answers.Rdata')
  newAnswers <- makeTestObjects()
  diffs <- abs(answers - newAnswers)
  max_diff_ind <- which(diffs == max(diffs) )[1]
  cat('Max Difference = ', max(diffs)[1], '\nMeasure = ', names(newAnswers)[max_diff_ind], '\n')
}