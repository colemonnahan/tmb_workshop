## Lecture 5
library(TMB)
library(ggplot2)
library(INLA)
library(reshape2)
## Load data and prepare for model
load('tmb_models/EBS_pollock_data.rda')
pollock <- EBS_pollock_data
pollock <- subset(pollock, year==2014)
## Matrix of pairwise differences
dd <- as.matrix(dist(pollock[,c('lat', 'long')]))
n <- length(pollock$catch)
pollock$lat <- pollock$lat-mean(pollock$lat)
pollock$long <- pollock$long-mean(pollock$long)

ggplot(pollock, aes(long, lat, size=(catch), color=catch==0)) + geom_jitter(alpha=.5)
compile("tmb_models/cpue_spatial.cpp")
dyn.load(dynlib("tmb_models/cpue_spatial"))
pars = list(intercept=3, theta=-5, beta_lat=0, beta_lon=0,
            logsigma=0, logsigma_space=1, a=1, u=rep(0, n))
data = list(likelihood=2, y=pollock$catch, lat=pollock$lat,
            lon=pollock$long, dd=dd)
obj = MakeADFun(data=data, parameters=pars, random='u', DLL="cpue_spatial")
#obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
opt <- with(obj, nlminb(opt$par, fn, gr)) # restart optimizer
rep1 <- obj$report()

### The SPDE approach to simplify the model. Using RINLA tools to create
### inputs to TMB. Based of Jim Thorson's function.
create.spde <- function(dat, n_knots, make.plot=FALSE, jitter=.3){
  loc_xy <- data.frame(x=dat$long, y=dat$lat)
  knots <- kmeans( x=loc_xy, centers=n_knots )
  loc_centers <- knots$centers
  ## loc_xy <- cbind(loc_xy, cluster=knots$cluster)
  ## ggplot(loc_xy, aes(x,y, col=factor(cluster))) + geom_point(size=.5)
  mesh <- inla.mesh.create( loc_centers, refine=TRUE)
  spde <- inla.spde2.matern( mesh )
  if(make.plot){
    png(paste0('mesh_', n_knots, '.png'), width=7, height=4, units='in', res=500)
    par(mar=.5*c(1,1,1,1))
    plot(mesh, main=NA, edge.color=gray(.7))
    points( jitter(dat$longitude, amount=jitter), jitter(dat$latitude, amount=jitter), cex=1, pch='.', col=rgb(0,0,0,.3))
    points( loc_centers, cex=.5, pch=20, col="red")
    dev.off()
  }
  return(list(mesh=mesh, spde=spde, cluster=knots$cluster, loc_centers=loc_centers))
}
spde <- create.spde(dat=pollock, n_knots=50, make.plot=FALSE)

## Make spde inputs for TMB
data <- list(likelihood=2, y=pollock$catch, lat=pollock$lat,
             lon=pollock$long, site=spde$cluster,
             M0=spde$spde$param.inla$M0, M1=spde$spde$param.inla$M1,
             M2=spde$spde$param.inla$M2)
pars <- list(intercept=3, theta=-5, beta_lat=0, beta_lon=0, logsigma=0,
             logsigma_space=1, u=rnorm(spde$mesh$n,0,1), logkappa=-3)
compile("tmb_models/cpue_spatial_spde.cpp")
dyn.load(dynlib("tmb_models/cpue_spatial_spde"))

## We can use map to turn off spatial component
map <- list(logkappa=factor(NA), logsigma_space=factor(NA), u=factor(rep(NA, length(pars$u))))
obj = MakeADFun(data=data, parameters=pars, 
                map=map, DLL="cpue_spatial_spde")
#obj$env$beSilent()
opt2 <- with(obj, nlminb(par, fn, gr))
opt2 <- with(obj, nlminb(opt2$par, fn, gr)) # restart
rep2 <- obj$report()

## Rerun with spatial SPDE approach
obj = MakeADFun(data=data, parameters=pars, random='u', 
                 DLL="cpue_spatial_spde")
#obj$env$beSilent()
opt3 <- with(obj, nlminb(par, fn, gr))
opt3 <- with(obj, nlminb(opt3$par, fn, gr)) # restart
opt3 <- with(obj, nlminb(opt3$par, fn, gr)) # restart
rep3 <- obj$report()
rep3$u

library(reshape2)
## Get residuals for all three models
pollock$resids.spde <- (rep3$pred-log(data$y))/rep3$sigma
pollock$resids.space <- (rep1$pred-log(data$y))/rep1$sigma
pollock$resids.nospace <- (rep2$pred-log(data$y))/rep2$sigma
resids <- melt(pollock, id.vars = c("lat", "long"), measure.vars = c("resids.space", "resids.nospace", "resids.spde"))
ggplot(resids, aes(long, lat, size=abs(value), color=value)) +
  geom_point(alpha=.5)+ scale_color_gradient(low='red', high='black') +
  facet_wrap('variable') + scale_size_area(max_size=4)



### Part 2: Bayesian inference
## We return to the swallows mark-recapture model and refit as a Bayesian
## model using Stan
library(TMB)
data <- readRDS('tmb_models/swallows.RDS')
pars <- list(sigmayearphi=0,
         sigmaphi=0,
         sigmap=0,
         a=rep(3, data$K-1),
         a1=1,
         b0=rep(1,4),
         b1=rep(0,4),
         fameffphi_raw=rep(0,data$nfam),
         fameffp_raw=rep(0,data$nfam),
         yeareffphi_raw=rep(0,4))

## Note this model is setup to be run as a Bayesian model. We can ignore
## that for now
compile("tmb_models/swallows.cpp")
dyn.load(dynlib("tmb_models/swallows"))
obj <- MakeADFun(data, pars, random=c("fameffphi_raw", "fameffp_raw", "yeareffphi_raw"), DLL='swallows')
opt <- with(obj, nlminb(par, fn, gr))
opt <- with(obj, nlminb(opt$par, fn, gr)) # restart

library(tmbstan)
library(shinystan)

## This plugs our TMB model into Stan and returns a fitted MCMC object,
## just like if you had a model in Stan
fit <- tmbstan(obj, chains=1, iter=1000) # short chain for demo
fit

launch_shinystan(fit)

## ## This takes a long time to run
## fit1 <- tmbstan(obj, laplace=FALSE)
## fit2 <- tmbstan(obj, laplace=TRUE)
