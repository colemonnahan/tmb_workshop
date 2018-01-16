# Exercise: Asymptotic uncertainty
compile("tmb_models/growth_example.cpp")
dyn.load(dynlib("tmb_models/growth_example"))
dat <- read.table("tmb_models/bevholt.dat", header=TRUE)
#data <- list(SSB=dat$ssb,logR=dat$logR)
data <- readRDS('tmb_models/growth_data.RDS')
parameters <- list(logsigma=0, logLinf=3, logk=-2)
#parameters <- list(logA=0, logB=0, logsigma=0)

obj <- MakeADFun(data,parameters, DLL="growth_example")

opt <- with(obj, nlminb(par, fn, gr))
opt
