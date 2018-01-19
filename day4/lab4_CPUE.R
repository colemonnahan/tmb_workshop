library(plyr)
library(TMB)
library(ggplot2)

### Exercise 1: CPUE standardization of pollock
## Catch rate standardization (catch per unit effort) for Alaska
## Pollock. Modified with permission from lectures at
## https://github.com/James-Thorson/2016_Spatio-temporal_models

## Data taken on 2017/12/20 from
## https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/data/EBS_pollock_data.rda


## Load data and prepare for model
load('tmb_models/EBS_pollock_data.rda')
pollock <- EBS_pollock_data
pollock <- subset(pollock, year==2014)
n <- length(pollock$catch)
ggplot(pollock, aes(long, lat, size=(catch), color=catch==0)) + geom_jitter(alpha=.5)
compile("tmb_models/cpue1.cpp" )
dyn.load( dynlib("tmb_models/cpue1") )
pars = list("intercept"=5, "theta"=0, "beta_lat"=0, "beta_lon"=0, "logsigma"=0)

## First with invgauss likelihood
data = list(likelihood=1, y=pollock$catch, lat=pollock$lat, lon=pollock$long)
obj = MakeADFun(data=data, parameters=pars, DLL="cpue1")
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
AIC1 <- 2*opt$objective-2*length(opt$par)

## Repeat with lognormal
data = list(likelihood=2, y=pollock$catch, lat=pollock$lat, lon=pollock$long)
obj = MakeADFun(data=data, parameters=pars, DLL="cpue1")
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
AIC2 <- 2*opt$objective-2*length(opt$par)

## Clear support for lognormal likelihood
c(AIC1, AIC2)

## Plot residuals
rep <- obj$report()
y <- data$y
yhat <- rep$pred
pollock$resids.nospace <- (yhat-log(y))/rep$sigma
ggplot(pollock, aes(long, lat, size=abs(resids), color=resids)) +
  geom_point(alpha=.5)+ scale_color_gradient(low='red', high='blue') +
  scale_size_area(max_size=4)

### Part 3, use cross validation to test lat/lon covariates
## First with pars turned off
compile("tmb_models/cpue2.cpp" )
dyn.load( dynlib("tmb_models/cpue2") )

K <- 10
partition <- sample( x=1:K, size=length(data$y), replace=TRUE )

data = list(likelihood=2, y=pollock$catch, lat=pollock$lat, lon=pollock$long)
predNLL <- rep(NA, K)
## Try different maps options. Which has the lowest mean predicted error? 
#map <- list(beta_lat=factor(NA), beta_lon=factor(NA))
map <- list()
for(k in 1:K){
  ## Resample which data points to exclude from NLL calculation (10%)
  data$out_sample <- ifelse(partition==k,0,1)
  ## Remake obj and solve
  obj = MakeADFun(data=data, parameters=pars, DLL="cpue2", map=map)
  obj$env$beSilent()
  opt <- with(obj, nlminb(par, fn, gr))
  ## The total NLL/total points = average NLL per point
  predNLL[k] <- obj$report()$pred_jnll/sum(partition==k)
}
mean(predNLL)
