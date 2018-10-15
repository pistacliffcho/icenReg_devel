# TESTING BOOTSTRAP CODE

library(icenReg)

simData = simIC_cluster()

fit_np = ic_np(simData[,c("l", "u")]  )
fit_par = ic_par(cbind(l, u) ~ x1 + x2, data = simData)
fit_bayes = ic_bayes(cbind(l,u) ~ x1 + x2, data = simData)
fit_sp = ic_sp(cbind(l,u) ~ x1 + x2, data = simData)

fit_np
fit_bayes
fit_par
fit_sp

ir_clustBoot(fit_bayes, simData$ID)
ir_clustBoot(fit_np, simData$ID)
ir_clustBoot(fit_par, simData$ID)
ir_clustBoot(fit_sp, simData$ID)

fit_par
fit_sp
