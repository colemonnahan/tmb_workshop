library(plyr)
library(TMB)
library(ggplot2)
###  von Bertalanffy growth of (simulated) fish
sample.lengths <- function(Nfish, n.ages, logLinf.mean=log(50), logLinf.sigma=.1,
                           logk.mean=log(.1), logk.sigma=.1, sigma.obs=.1,
                           t0=5, Ntime=40){
  sample.ages <- function(n.ages, t0, Ntime) {sample((t0+1):Ntime, size=n.ages, replace=FALSE)}
  sample.vbgf <- function(ages, Linf, k,  t0, sigma.obs){
    lengths <- Linf*(1-exp(-k*(ages-t0)))
    loglengths <- log(lengths)+ rnorm(n=length(lengths), mean=0, sd=sigma.obs)
    data.frame(ages=ages, loglengths=loglengths)
  }
    Linf.vec <- exp(logLinf.mean + rnorm(n=Nfish, 0, sd=logLinf.sigma))
    k.vec <- exp(logk.mean +rnorm(n=Nfish, mean=0, sd=logk.sigma))
    dat <- ldply(1:Nfish, function(i)
        cbind(fish=i, sample.vbgf(ages=sample.ages(n.ages, t0=t0, Ntime=Ntime),
              Linf=Linf.vec[i], k=k.vec[i], sigma.obs=sigma.obs, t0=t0)))
   return( dat)
}

## Generate data
N <- 20
seed <- 313415
set.seed(seed)
dat <- sample.lengths(Nfish=N, n.ages=5)
g <- ggplot(dat, aes(ages, exp(loglengths), group=fish)) +
    geom_point(alpha=.5, size=.1) + geom_line(alpha=.5) + xlab("Time (t)") + ylab("Length")
g + theme_bw()
data <- list(Nfish=N, Nobs=nrow(dat), loglengths=dat$loglengths,
             fish=dat$fish, ages=dat$ages)
saveRDS(file='tmb_models/growth_data.RDS', data)


## Exercise 1: Growth model with common parameters
data <- readRDS('tmb_models/growth_data.RDS')
parameters <- list(logsigma=0, logLinf=3, logk=-2)
data$ages_pred <- 6:40
compile('tmb_models/growth.cpp')
dyn.load(dynlib('tmb_models/growth'))
obj <- MakeADFun(data=data, parameters=parameters, DLL='growth')
obj$env$beSilent()
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
(mles <- rep$value)
ses <- rep$sd
plot(data$ages, exp(data$loglengths), cex=.5)
lines(data$ages_pred, mles, lwd=2)
lines(data$ages_pred, mles-1.96*ses, lty=2)
lines(data$ages_pred, mles+1.96*ses, lty=2)


## Exercise 2: Growth model with independent parameters
N <- 20
data <- readRDS('tmb_models/growth_data.RDS')
data$ages_pred <- 6:40
parameters <- list(logsigma=0, logLinf=rep(3, N), logk=rep(-2, N))
compile('tmb_models/growth2.cpp')
dyn.load(dynlib('tmb_models/growth2'))
obj <- MakeADFun(data=data, parameters=parameters, DLL='growth2')
obj$env$beSilent()
opt <- nlminb(obj$par, obj$fn, obj$gr)
max(abs(obj$gr(opt$par)))               # gradient not small
opt$message                             # needs more iterations
## Rerun starting from opt
opt <- nlminb(opt$par, obj$fn, obj$gr)
max(abs(obj$gr(opt$par)))               # gradient small

## Get the MLEs and SEs from the model.  Drop the 
## first parameter since it is logsigma
rep <- sdreport(obj)
ses <- sqrt(diag(rep$cov.fixed))[-1]
mles <- opt$par[-1]
## Plot them, one for each fish
df <- data.frame(fish=1:N, par=rep(c('logLinf', 'logK'), each=N),
                 mles=mles, lwr=mles-1.96*ses,
                 upr=mles+1.96*ses)
ggplot(df, aes(x=fish, y=mles, ymin=lwr, ymax=upr)) + geom_pointrange() +
                 facet_wrap('par', scales='free', ncol=1)

