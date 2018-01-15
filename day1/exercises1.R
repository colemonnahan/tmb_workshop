library(TMB)

## Exercise 1: using TMB derivative
compile("TMB_models/polynomial.cpp")
dyn.load(dynlib("TMB_models/polynomial"))
obj <- MakeADFun(data=list(), parameters=list(x=0), DLL='polynomial')
obj$fn(5)
obj$gr(5)

## Solution
x <- seq(-3,2, len=200)
y <- g <- rep(NA, 200)
 for(i in 1:length(x)) {
  y[i] <- obj$fn(x[i])
  g[i] <- obj$gr(x[i])
 }
png("polynomial.png", width=5.5, height=2.5, units='in', res=500)
par(mfrow=c(1,2), tck=-.02, mgp=c(1.5,.25,0), mar=c(3,3,1,.1))
plot(x,y, type='l', lwd=2, main='Height', ylab='fn')
plot(x,g, type='l', lwd=2, main='Derivative', ylab='gr')
abline(h=0)
dev.off()

## Exercise 2 solution: Poisson likelihood
k <- 4
loglike <- function(lambda) k*log(lambda)-lambda-log(factorial(k))
loglike(5.5)
dpois(x=4, lambda=5.5, log=TRUE)
## lgamma is better because it is more stable numerically
lgamma(k+1)
log(factorial(k))
ll <- seq(0.1, 15, len=1000)
plot(ll, -loglike(ll), type='l', lwd=2)



