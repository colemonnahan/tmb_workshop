library(TMB)
library(plyr)
library(ggplot2)

## Exercise 1: Poisson count data
## Simulate Poisson example
inputs <- function(mu, sigma, nsites, nreps){
  D <- rnorm(nsites, mu, sigma)
  sites <- rep(1:nsites, each=nreps)
  C <- rpois(nsites*nreps, D[sites])
  data <- list(C=C, sites=sites)
  pars <- list(D=rep(mu,nsites), mu=mu, logsigma=log(sigma))
  return(list(data=data, pars=pars))
}

compile("tmb_models/poisson.cpp")
dyn.load(dynlib("tmb_models/poisson"))

## The original one with little information
set.seed(3432)
x <- inputs(20, 5, 4, 2)
obj <- MakeADFun(data=x$data, parameters=x$pars, random="D", DLL='poisson')
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
## Can sometimes be -Inf for lower, so constrain it
prof1 <- tmbprofile(obj, 'logsigma', parm.range=c(0,5))
confint(prof1)
par(mfrow=c(1,3))
plot(prof1)

## Add more sites and replicates
set.seed(3432)
x <- inputs(20, 5, 50, 5)
obj <- MakeADFun(data=x$data, parameters=x$pars, random="D", DLL='poisson')
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
prof2 <- tmbprofile(obj, 'logsigma')
plot(prof2)

## Even more data
set.seed(3432)
x <- inputs(20, 5, 5000, 50)
obj <- MakeADFun(data=x$data, parameters=x$pars, random="D", DLL='poisson')
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
prof3 <- tmbprofile(obj, 'logsigma')
plot(prof3)

## Make table
table1 <- exp(rbind(confint(prof1), confint(prof2), confint(prof3)))
dimnames(table1)[[1]] <- paste('scenario', 1:3)
table1




### Model 2: Hierarchical growth with random effects
### Simulation growth functions
N <- 20
data <- readRDS('tmb_models/growth_data.RDS')
data$ages_pred <- 6:40

## Refit model with no pooling
parameters <- list(logsigma=-2.5, logLinf=rep(3, N), logk=rep(-2, N))
compile('tmb_models/growth2.cpp')
dyn.load(dynlib('tmb_models/growth2'))
obj <- MakeADFun(data=data, parameters=parameters, DLL='growth2')
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt <- nlminb(opt$par, obj$fn, obj$gr) # restart
rep <- sdreport(obj)
ses <- sqrt(diag(rep$cov.fixed))[-1]
mles <- opt$par[-1]
## Plot them, one for each fish
df <- data.frame(fish=1:N, par=rep(c('logLinf', 'logK'), each=N),
                 mles=mles, lwr=mles-1.96*ses,
                 upr=mles+1.96*ses)

## Fit random effect model ( partial pooling)
parameters <- list(logsigma_obs=-2.5, logLinf=rep(3, N), 
                   logk=rep(-2, N),
                   logsigma_logLinf=0, mean_logLinf=3, 
                   logsigma_logk=0, mean_logk=-2)
compile('tmb_models/growth3.cpp')
dyn.load(dynlib('tmb_models/growth3'))
obj <- MakeADFun(data=data, parameters=parameters,
                 random=c('logk', 'logLinf'), DLL='growth3')
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
rep2 <- sdreport(obj)

## Make a quick plot of comparisons
df2 <- data.frame(fish=1:N, par=rep(c('logLinf', 'logK'), each=N),
                 mles=rep2$value, lwr=rep2$value-1.96*rep2$sd,
                 upr=rep2$value+1.96*rep2$sd)
df.all <- rbind(cbind(model='no pooling', df), cbind(model='partial pooling', df2))
ggplot(df.all, aes(x=fish, y=mles, ymin=lwr, ymax=upr, color=model)) + geom_point() +
                  facet_wrap('par', scales='free', ncol=1)
