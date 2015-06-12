
library(icenReg)
library(intcox)
library(straweib)

set.seed(1)

#Processing functions
#Processes the llkArray by subtracting off the max likelihood for each simulated dataset
process_llk_array <- function(llkArray){
	dims <- dim(llkArray)
	I <- dims[2]
	J <- dims[3]
	
	outArray <- llkArray
	for(i in 1:I){
		for(j in 1:J){
			maxllk <- max(llkArray[,i,j])
			outArray[,i,j] <- outArray[,i,j] - maxllk
		}
	}
	return(outArray)
}

getRowMeans <- function(sumArray){
	#assumes 3 algorithms, 4 sample sizes
	ans <- matrix(nrow = 3, ncol = 4)
	for(i in 1:4){
		ans[1,i] <- mean(sumArray[1,,i])
		ans[2,i] <- mean(sumArray[2,,i])
		ans[3,i] <- mean(sumArray[3,,i])
	}
	rownames(ans) <- c('icph', 'intcox_default', 'intcox_strict')
	
	return(ans)
}



###				ACTUAL CODE FOR RUNNING SIMULATIONS/EXAMPLES

#		Code for demonstrating problem with intcox in tutorial
data(tooth24)
tooth24$rightNA <- tooth24$right				#Preprocessing mentioned in tutorial
tooth24$rightNA[tooth24$right == Inf] <- NA		
intcox_fit <- intcox(Surv(left, rightNA, type = 'interval2') ~ sex + dmf, data = tooth24)

icenReg_fit <- ic_sp(Surv(left, right, type = 'interval2') ~ sex + dmf, data = tooth24, bs_samples = 0)

intcox_fit$loglik
icenReg_fit$final_llk
#Very different final likelihoods!!

intcox_fit$coef
icenReg_fit$coef
#Different final regression coefficients


# Code for table of algorithm speed (ic_ph)

# Note: this will probably take about 10 hours! For a snap shot, set numClocks = 4

numClocks = 100					#Number of simulations per sample size
ns = c(100, 500, 1000, 5000)	#Sample size


timeArray <- array(dim = c(3,  numClocks, length(ns)), 
						 dimnames = list(algorithm = c('ic_ph', 'intcox_st', 'intcox_df'), 
						 values = NULL, sampleSizes = ns) )
llkArray <- timeArray
# The above is used to store the outcome

st_eps = 10^-6
# stricter stopping criterion used for intcox


#####			RUNNING AND PROCESSING OF TIMING OF ALGORITHMS
for(i in 1:length(ns)){
	for(j in 1:numClocks){
		cat('\nstarting ij = ', i, j, '\n')
#		Uncomment the above if you want to see how the simulation is going
		simdata <- simIC_weib(n = ns[i])
		is.inf <- simdata$u == Inf
		simdata$u[is.inf] <- NA		#intcox does not accept Inf, but works with NA
		timeArray[1, j, i] <- system.time(myFit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata, bs_samples = 0))[1]
		timeArray[2, j, i] <- system.time(intcoxFit <- intcox(Surv(l, u, type = 'interval2') ~ x1 + x2 , data = simdata))[1]
		timeArray[3, j, i] <- system.time(intcox_stFit <- intcox(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata, epsilon = st_eps, itermax = 100000))[1]
		llkArray[1, j ,i] <- myFit$final_llk
		llkArray[2, j, i] <- intcoxFit$loglik
		llkArray[3, j, i] <- intcox_stFit$loglik
	}
}


disc_timeMatrix <- matrix(nrow = numClocks, ncol = length(ns))
colnames(disc_timeMatrix) <- ns
# Used to save the outcomes with a discrete censoring
# Only need a matrix because intcox fails too often with discrete censoring, so only need to review one algorithm

icpo_timeMatrix <- disc_timeMatrix
disc_icpo_timeMatrix <- disc_timeMatrix
#These are for saving the times of the ic_po algorithm (no other algorithm to compare with)


for(i in 1:length(ns)){
	for(j in 1:numClocks){
		cat('ij = ', i, j, '\n')
		simdata <- simIC_weib(n = ns[i], rndDigits = 1)
		disc_timeMatrix[j,i] <- system.time(ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata, bs_samples = 0))[1]
	}
}

for(i in 1:length(ns)){
	for(j in 1:numClocks){
		cat('ij = ', i, j, '\n')
		simdata <- simIC_weib(n = ns[i], model = 'po')
		icpo_timeMatrix[j,i] <- system.time(ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata, bs_samples = 0, model = 'po'))[1]
	}
}

for(i in 1:length(ns)){
	for(j in 1:numClocks){
		simdata <- simIC_weib(n = ns[i], rndDigits = 1, model = 'po')
		disc_icpo_timeMatrix[j,i] <- system.time(ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata, bs_samples = 0, model = 'po'))[1]
	}
}


######			DISPLAY OF COMPARISONS OF TIMING OF ALGORITHMS
diffArray <- process_llk_array(llkArray)
getRowMeans(diffArray)		#Average difference between MLE and final log likelihood (icenReg:ph + intcox)
getRowMeans(timeArray)		#Average time (icenReg:ph + intcox)

colMeans(disc_timeMatrix)	#Average time for discrete censoring (icenReg:ph)

colMeans(icpo_timeMatrix)		#Average time (icenReg:po)

colMeans(disc_icpo_timeMatrix)		#Average time for discrete censoring (icenReg:po)


# code for "Example Analysis section"

library(doParallel)
myCluster <- makeCluster(3, type = 'FORK')
registerDoParallel(myCluster)

data(mdata)
sp_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ grp, 
	data = mdata, bs_samples = 500, useMCores = T)
summary(sp_fit)

newdata <- data.frame(grp = c('ce', 'ge'))
rownames(newdata) = c('Conventional Environment', 'Germ Free Environment')
plot(sp_fit, newdata)

sCurves <- getSCurves(sp_fit, newdata)





# clocking ic_par against survreg
K = 100
icenRegTime <- numeric()
survRegTime <- numeric()
lk_df <- numeric()
n <- 100000
for(i in 1:K){
	simdata <- simIC_weib(n = n, inspections = 20, rndDigits = 2)
	simdata$l <- simdata$l + 0.001		
	#survreg cannot handle left side = 0, even for general interval censoring
	is.inf <- simdata$u == Inf
	simdata$u[is.inf] <- 10^8
	#neither can it handle  Inf
	
	icenRegTime[i] <- system.time(icenRegFit <- ic_par(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata) )[1]
	survRegTime[i] <- system.time( survRegFit <- survreg(Surv(l, u, type = 'interval2') ~ x1 + x2, data = simdata) )[1]
	lk_df[i] <- icenRegFit$final_llk - survRegFit$loglik[2]
	if( (icenRegFit$final_llk - survRegFit$loglik[2]) < -1){cat("weird data found"); weird_data <- simdata}	
	if(icenRegFit$final_llk == -Inf) {bad_data <- simdata; stop('bad data found')}
}



#test