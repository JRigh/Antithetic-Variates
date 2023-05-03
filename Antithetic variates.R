#------------------------------------
# Antithetic variates for Monte Carlo
# variance reduction in R
#------------------------------------

# we want to estimate int from 0 to 1 of exp( x^2 ) dx

# plot the function
g = function(x) {exp( x^2)} # function g(x)
x <- seq(from = 0, to = 1, by = 0.01)
plot(x, g(x), type = "l", ylab = "g(x)", xlab = "x", col = 'darkred', lwd = 2,
     main = "g(x) = exp( x^2 ) on the interval [0,1]")

integrate(g, lower = 0, upper = 1)
# 1.462652 with absolute error < 1.6e-14

m <- 10000 # number of simulations
set.seed(1986) # for replication

# MC crude estimation
u <- runif(m)
gx_c <- exp(u^2)
I_hat_crude <- mean(gx_c)
I_hat_crude_2 <- (1/m) * sum(gx_c)
se_I_hat_crude <- sqrt(var(gx_c) / m)
se_I_hat_crude_2 <- sqrt(((1/m)^2)*sum((gx_c - I_hat_crude_2)^2 ))

# Antithetic estimation
u <- runif(m/2)
gx1 <- exp(u^2)
gx2 <- exp((1-u)^2)
gx <-(gx1+gx2)/2  # length(gx) [1] 5000
I_hat_anti <- mean(gx)
I_hat_anti_2 <- 2*(1/m) * sum(gx)
se_I_hat_anti <- sqrt(var(gx)/(m/2) )
se_I_hat_anti_2 <- sqrt((4* (1/m)^2)*sum((gx - I_hat_anti_2)^2 ))

# results  
results <- matrix(c(I_hat_crude,I_hat_anti,se_I_hat_crude, se_I_hat_anti),nrow=2)
colnames(results) <- c("estimate", "sd")
rownames(results) <- c("crude MC", "Antithetic")
results
#            estimate          sd
# crude MC   1.463267 0.004767316
# Antithetic 1.463597 0.002362512

# plot convergence of MC estimators
seq = 1:(m/2)
plot(cumsum(gx)/seq, type="l", ylab="MC estimators", col="firebrick3",
     xlab = "m/2", main = "Convergence of the MC estimators", lwd = 2)
lines(cumsum(gx_c[1:(m/2)])/seq, type="l", col="black",lwd = 2)
legend(3000, 1.52, legend=c("Crude MC estimator (black)", "Antithetic estimator (red)"),
       col=c('black', 'firebrick3'), lty=1, cex=0.7)

#----
# end
#---