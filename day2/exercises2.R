## Exercise: Add logsigma as parameter
library(TMB)
compile("tmb_models/bevholt2.cpp")
dyn.load(dynlib("tmb_models/bevholt2"))
dat <- read.table("tmb_models/bevholt.dat", header=TRUE)
data <- list(SSB=dat$ssb,logR=dat$logR)
parameters <- list(logA=0, logB=0, logsigma=0)
obj <- MakeADFun(data,parameters,DLL="bevholt2")
obj$env$beSilent() # silences console output
opt <- with(obj, nlminb(par, fn, gr))
rep <- sdreport(obj)
summary(rep, 'fixed')
summary(rep, 'report')


## Exercise: Asymptotic uncertainty
compile("tmb_models/bevholt3.cpp")
dyn.load(dynlib("tmb_models/bevholt3"))
dat <- read.table("tmb_models/bevholt.dat", header=TRUE)
data <- list(SSB=dat$ssb,logR=dat$logR, 
             SSB_pred=seq(20000, 300000, len=1000))
parameters <- list(logA=0, logB=0, logsigma=0)

obj <- MakeADFun(data,parameters, DLL="bevholt3")
obj$env$beSilent() # silences console output
opt <- with(obj, nlminb(par, fn, gr))
rep <- sdreport(obj) # AFTER optimizing model

## SEs of parameters. First, directly from the Hessian
sqrt(diag(solve(obj$he(opt$par))))
## Or can use this command
summary(sdreport(obj), 'fixed')

## SEs of predictions
dev.off()
ssb.df <- data.frame(SSB=data$SSB_pred, logR=rep$value, se=rep$sd)
plot(data$SSB, data$logR, xlab='SSB', ylab='logR')
lines(data$SSB_pred, ssb.df$logR, type='l', lwd=2)
lines(data$SSB_pred, ssb.df$logR+1.96*ssb.df$se, type='l')
lines(data$SSB_pred, ssb.df$logR-1.96*ssb.df$se, type='l')

## Exercise: Bootstrapping on predicted values
Nboot <- 200
boot.results <- matrix(NA, nrow=Nboot, ncol=2+1000)
N <- length(data$SSB)
for(i in 1:Nboot){
  ind <- sample(1:N, size=N, replace=TRUE)
  data2 <- list(SSB=data$SSB[ind], logR=data$logR[ind],
                SSB_pred=seq(20000, 300000, len=1000))
  obj.temp <- MakeADFun(data2,parameters,DLL="bevholt3")
  obj.temp$env$beSilent() # silences console output
  opt.temp <- with(obj.temp, nlminb(par, fn, gr))
  boot.results[i,] <- sdreport(obj.temp)$value
}
## Plot them
par(mfrow=c(1,1))
plot(data$SSB, data$logR, pch=16, col=1, ylab='logR', xlab='SSB')
apply(boot.results, 1, function(x)
  lines(data$SSB_pred, x[-(1:2)], col=rgb(0,0,0,.1)))
lines(data$SSB_pred, ssb.df$logR+1.96*ssb.df$se, type='l', col=2, lwd=2)
lines(data$SSB_pred, ssb.df$logR-1.96*ssb.df$se, type='l', col=2, lwd=2)
lines(data$SSB_pred, ssb.df$logR, type='l')
points(data$SSB, data$logR, pch=16, col=1)
