library(TMB)
library(ggplot2)
## Catch rate standardization (catch per unit effort) for Alaska
## Pollock. Modified with permission from lectures at
## https://github.com/James-Thorson/2016_Spatio-temporal_models

## Data taken on 2017/12/20 from
## https://github.com/nwfsc-assess/geostatistical_delta-GLMM/blob/master/data/EBS_pollock_data.rda

setwd('../TMB_models')

## Load data and prepare for model
load('EBS_pollock_data.rda')
pollock <- EBS_pollock_data
pollock <- subset(pollock, year==2014)
ggplot(pollock, aes(long, lat, size=(catch), color=catch==0)) + geom_jitter(alpha=.5)
compile( "cpue.cpp" )
dyn.load( dynlib("cpue") )

pars = list("intercept"=0, "theta"=0, "logsigma"=0)
data = list(likelihood=10, y=pollock$catch)
obj = MakeADFun(data=data, parameters=pars, DLL="cpue")
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))

## TODO: fit AIC for each likelihood and show the cross validation
## technique

for(j in 1:3){
  Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
  Data = list( "Options_vec"=j-1, "y"=CPUE, "X_ij"=X, predTF_i=ifelse(Partition_i==k,1,0) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="cpue")
  Obj$env$beSilent()

  # Step 3 -- test and optimize
  Opt = nlminb( start=Obj$env$last.par.best, objective=Obj$fn, gradient=Obj$gr )

  # Extract stuff
  Report = Obj$report()
  PredNLL_kj[k,j] = Report$pred_jnll / sum(Partition_i==k)
  if( PredNLL_kj[k,j]>10 ) stop("Check results")
}


