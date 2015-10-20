library(icenReg)
n = 500
simdata <- simIC_weib(n, model = 'ph', b1 = -.1, b2 = .5, prob_cen = 0)

fit <- ic_sp(cbind(l, u) ~ x1 + x2, maxIter = 100, data = simdata)
#fit <- ic_sp(cbind(l, u) ~ x1 + x2, maxIter = 500, data = simdata)
#fit <- ic_sp(cbind(l, u) ~ x1 + x2, maxIter = 500, model = 'po', data = simdata)

iterRange <- 5:35
fits <- list()
baseUpdates <- 15
colChanges <- 5
useGA = TRUE
for(i in seq_along(iterRange)){
  thisIter <- iterRange[i]
  fits[[i]] <- ic_sp(cbind(l, u) ~ x1 + x2, data = simdata, 
                     maxIter = thisIter, baselineUpdates = baseUpdates,
                     useGA = useGA)
}

p_diffs = list()
minVal = Inf
maxVal = -Inf
for(i in 2:length(fits)){
  p_diffs[[i-1]] <- fits[[i]]$p_hat - fits[[i-1]]$p_hat
  minVal = min(c(p_diffs[[i-1]], minVal) )
  maxVal= max(c(p_diffs[[i-1]], maxVal) )
}

plot(1:length(p_diffs[[1]]), p_diffs[[1]], type = 'l', ylim = c(minVal, maxVal),
     xlab = "parameter", ylab = "Change in parameter from last iteration", main = "Converge-O-Gram")
for(i in 2:length(p_diffs)){
  thisCol = ceiling(i/colChanges)
  lines(1:length(p_diffs[[1]]), p_diffs[[i]], col = thisCol)
}

numcols = 1:(ceiling(length(fits) / colChanges))
legend('topleft', legend = numcols, col = numcols, lty = 1)

lines(c(0, 1000000), c(0,0), lwd = 2.5)

last_diffs <- numeric()
for(i in 1:length(p_diffs)){
  last_diffs[i] <- p_diffs[[i]][n]
}

last2_diffs <- numeric()
for(i in 1:length(p_diffs)){
  last2_diffs[i] <- p_diffs[[i]][n-1]
}