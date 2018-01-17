## Lecture 2.3 code
library(TMB)
library(ggplot2)

## Simulate Poisson example
## set.seed(3432)
## mu <- 20
## sigma <- 5
## D <- rnorm(4, mu, sigma)
## sites <- rep(1:4, each=2)
## C <- rpois(8, D[sites])
C <- c(10, 11, 24, 21, 16, 13, 14, 10)
sites <- c(1, 1, 2, 2, 3, 3, 4, 4)

compile("tmb_models/poisson.cpp")
dyn.load(dynlib("tmb_models/poisson"))
data <- list(C=C, sites=sites)
pars <- list(D=rep(20,4), mu=20, logsigma=1)
obj <- MakeADFun(data=data, parameters=pars, DLL='poisson')
obj$fn()
opt <- with(obj, nlminb(par, fn, gr))
opt$par
prof <- tmbprofile(obj, 'logsigma')
plot(prof, main="Likelihood profile for logsigma")

## Look at NLL as variance goes to 0
-dnorm(0,0,1e-3,T)
-dnorm(0,0,1e-6,T)
-dnorm(0,0,1e-9,T)

## Remake with random effects turned on
obj <- MakeADFun(data=data, parameters=pars, random=c('b',"D"), DLL='poisson')
obj$fn()
obj$gr()
opt <- with(obj, nlminb(par, fn, gr))
prof <- tmbprofile(obj, 'logsigma', parm.range=c(-2,3))
plot(prof, main="Likelihood profile for logsigma")


