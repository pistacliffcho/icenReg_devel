library(icenReg)
library(ICsurv)
data("Hemophilia")
data(miceData)

hemoReady <- Hemophilia
hemoReady$R[ hemoReady$d3 == 1] <- Inf

cnts <- bayesControls(useMLE_start = F, 
                      acceptRate = 0.33,
                      updateChol = T, 
                      burnIn = 5000, 
                      iterationsPerUpdate = 100, 
                      thin = 1,
                      samples = 5000)

pFit <- ic_par(cbind(L, R) ~ Low + Medium + High, data = hemoReady)
bFit <- ic_bayes(cbind(L, R) ~ Low + Medium + High, data = hemoReady, controls = cnts)
summary(bFit) 
#plot(bFit$samples)
#plot(bFit$logPosteriorDensities, type = 'l')
plot(pFit, col = 'blue', lwd = 2)
lines(bFit, col = 'purple', lwd = 2)

plot(bFit, col = 'blue', lwd = 2)
lines(pFit, col = 'purple', lwd = 2)

newVals <- data.frame(Low = c(0, 1), Medium = c(1,0), High = c(0,0))
rownames(newVals) <- c('Medium', 'Low')

plot(pFit, newVals, lwd = 3, col = 'brown', plot_legend = F)
#lines(pFit, newVals, lwd = 2, col = 'red')


cnts1 <- bayesControls(samples = 10^4 * 5)
cnts2 <- bayesControls(samples = 10^4 * 5, useMLE_start = F, updateChol = T)
cnts3 <- bayesControls(samples = 10^4 * 5, updateChol = T)

fit1 <- ic_bayes(cbind(l, u) ~ grp, data = miceData, controls = cnts1,
                                dist = 'lnorm')
fit2 <- ic_bayes(cbind(l, u) ~ grp, data = miceData, controls = cnts2,
                                dist = 'lnorm')

mcmc1 <- mcmc(fit1$samples)
mcmc2 <- mcmc(fit2$samples)

summary(mcmc1)
summary(mcmc2)

effectiveSize(mcmc1)
effectiveSize(mcmc2)

quartz()
plot(mcmc1)
plot(mcmc2)


priorFxn <- function(betas){
  return(sum(dnorm(betas, log = T)))
}

fit3 <- ic_bayes(cbind(l, u) ~ grp, data = miceData, priorFxn = priorFxn, controls=  cnts2)
mcmc3 <- mcmc( fit3$samples )
plot(mcmc3)


data("IR_diabetes")
cnts1 <- bayesControls(samples = 10^4 * 5)
fit1 <- ic_bayes(cbind(left, right) ~ gender, data = IR_diabetes, controls = cnts1)
fit2 <- ic_bayes(cbind(left, right) ~ gender, data = IR_diabetes, controls = cnts2)
fit3 <- ic_bayes(cbind(left, right) ~ gender, data = IR_diabetes, controls = cnts3)

mcmc1 <- mcmc(fit1$samples)
mcmc2 <- mcmc(fit2$samples)
mcmc3 <- mcmc(fit3$samples)

summary(mcmc1)
summary(mcmc2)

effectiveSize(mcmc1)
effectiveSize(mcmc2)
effectiveSize(mcmc3)
quartz()
plot(mcmc1)
plot(mcmc2)
