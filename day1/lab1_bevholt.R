## Beverton-Hold stock recruit function. Modified from A. Nielson's
## lectures: http://www.nielsensweb.org/bergen/
library(TMB)

compile("tmb_models/bevholt.cpp")
dyn.load(dynlib("tmb_models/bevholt"))

dat <- read.table("tmb_models/bevholt.dat", header=TRUE)
data <- list(SSB=dat$ssb,logR=dat$logR)
parameters <- list(logA=0, logB=0)

obj <- MakeADFun(data,parameters,DLL="bevholt")
obj$env$beSilent() # silences console output
obj$fn()
obj$gr()

### --------------------------------------------------
### Exercise answers
opt <- nlminb(obj$par,obj$fn,obj$gr)
## The MLEs:
logA <- opt$par[1]
logB <- opt$par[2]
plot(dat$ssb,dat$logR)
ssb <- seq(0, 300000, len=1000)
ssb <- data$SSB
## Predicted logR:
logR <- logA+log(ssb)-log(1+exp(logB)*ssb)
lines(ssb, logR, lwd=3)

## Reoptimize model from a grid of starting points
x1 <- seq(-5, 10, len=100)
x2 <- seq(-15,0, len=100)
inits <- expand.grid(logA=x1,logB=x2)
nlls <- evals <- rep(NA, len=nrow(inits))
for(i in 1:nrow(inits)){
  opt.temp <- nlminb(inits[i,],obj$fn,obj$gr)
  evals[i] <- opt.temp$iterations
  nlls[i] <- opt.temp$objective
}
## These should all be the same:
summary(nlls)
summary(evals)

## Make a contour of the negative log-likelihood surface
o <- 6 # number of SEs away from MLE
logA.se <- sqrt(diag(sdreport(obj)$cov.fixed))[1]
logB.se <- sqrt(diag(sdreport(obj)$cov.fixed))[2]
x1 <- seq(logA-o*logA.se, logA+o*logA.se, len=100)
x2 <- seq(logB-o*logB.se, logB+o*logB.se, len=100)
z <- sapply(x1, function(logA) sapply(x2, function(logB)
  obj$fn(c(logA,logB))))
contour(x1,x2,z, nlevels=100)
points(logA, logB, pch=16, col=2)
covar <- sdreport(obj)$cov.fixed
## Add the confidence ellipse and standard errors
require(ellipse)
CI.region <- ellipse(x=covar, centre=c(logA, logB))
lines(CI.region, col=2)
lines(x=c(logA-1.96*logA.se, logA+1.96*logA.se), y=c(logB,logB), col=4)
lines(y=c(logB-1.96*logB.se, logB+1.96*logB.se), x=c(logA,logA), col=4)

## Simulation testing the BH model
Nsim <- 2000
mles <- ses <- coverage <- matrix(NA, nrow=Nsim, ncol=2)
logR.expected <- logA+log(data$SSB)-log(1+exp(logB)*data$SSB)
## Note: sigma=0.4 is fixed in the template and matches here
for(i in 1:Nsim){
  error <- rnorm(n=length(logR.expected), mean=0, sd=0.4)
  data$logR <- logR.expected+error
  obj <- MakeADFun(data,parameters,DLL="bevholt")
  obj$env$beSilent()
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  mles[i,] <- opt$par
  ses[i,] <- sqrt(diag(sdreport(obj)$cov.fixed))
  ## Coverage is whether the confidence interval covers the true
  ## parameter.
  coverage[i,1] <-
    mles[i,1]-1.96*ses[i,1] < logA &
    logA < mles[i,1]+1.96*ses[i,1]
  coverage[i,2] <-
    mles[i,2]-1.96*ses[i,2] < logB &
    logB < mles[i,2]+1.96*ses[i,2]
}

par(mfrow=c(1,2))
breaks <- 20
hist(mles[,1], main='logA', breaks=breaks)
abline(v=logA, col=2)
hist(mles[,2], main='logB', breaks=breaks)
abline(v=logB, col=2)

## Check coverage
100*apply(coverage, 2, mean)
