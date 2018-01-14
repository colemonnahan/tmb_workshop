library(TMB)

compile("TMB_models/polynomial.cpp")
dyn.load(dynlib("TMB_models/polynomial"))
obj <- MakeADFun(data=list(), parameters=list(x=0), DLL='polynomial')
obj$fn(5)
obj$gr(5)

plot_derv <- function(x, eps=.1, col=2){
  y <- obj$fn(x)
  g <- obj$gr(x)
  points(x,y,pch=16)
  xx <- c(x-eps, x+eps)
  yy <- c(y-eps*g, y+eps*g)
  lines(xx,yy, col=col, lwd=2.5)
}

for(i in 1:3){
  png(paste0("dervs",i,".png"), width=4.5, height=2.5, units='in', res=500)
  par(mfrow=c(1,1), mgp=c(1.5, .25, 0), tck=-.02, mar=c(3,3,.1,.1))
  plot(x,y, main=NA, type='n', ylab='y=f(x)')
  if(i==3) lines(x,y, lwd=2, col=gray(.5))
  points(1.13, obj$fn(1.13), pch=16)
  if(i>1) plot_derv(1.13)
  if(i==3){
    plot_derv(-2.5)
    plot_derv(0.13)
    plot_derv(-2)
    plot_derv(-1)
    opt <- nlminb(start=1, objective=obj$fn, gradient=obj$gr, lower=-3, upper=3)
    plot_derv(opt$par, col=4)
  }
  dev.off()
}


## Part 2: Fitting linear model: y=a+b*x
x <- c(1.87, 1.96, 1.39, 2.24, 2.33, 2.24, 2.67, 2.47, 1.35, 2.00)
y <- c(2.47, 2.42, 2.2, 2.72, 2.65, 2.5, 2.85, 2.77, 2.28, 2.45)
plot(x,y)

## Use the built-in lm function
fit1 <- lm(y~x)
abline(fit1, lwd=2)
coef(fit1)

## Manual NLL function
nll <- function(pars){
  ## Extract parameters from vector
  intercept <- pars[1]
  slope <- pars[2]
  sd <- exp(pars[3])
  ## Predict y given parameters
  mu <- intercept + slope*x
  ## Calculate log likelihood by hand
  nll <- -sum(dnorm(y, mu, sd, log=T))
  return(nll)
}

## Third way is to use TMB to calculate likelihood
compile("tmb_models/linmod.cpp")
dyn.load(dynlib("TMB_models/linmod"))
obj <- MakeADFun(data=list(x=x,y=y),
                 parameters=list(intercept=1, slope=0.5, logsd=-2),
                 DLL='linmod')

## Check that all three ways matches other ways
pars <- c(1.2, 0.6, -2)
nll(pars)
obj$fn(pars)
## TMB also provides gradient using automatic differentiation
obj$gr(pars)
fit.R <- nlminb(start=pars, objective=nll)
fit.tmb <- nlminb(start=pars, objective=obj$fn, gradient=obj$gr)
coef(fit1)
(mle <- fit.R$par)
## Should be close to zero:
obj$gr(fit.tmb$par)

## Show the improvement in NLL from intitial to MLE
obj$fn(pars)
obj$fn(mle)
plot(x,y)
abline(a=pars[1], b=pars[2], col='blue', lty=2, lwd=2)
abline(a=mle[1], b=mle[2], col='red', lty=3, lwd=2)

## What about the variance estimate?
summary(lm(y~x))
exp(fit.tmb$par[3])
## The sigma parameters do not match because lm() uses REML. We can
## recreate that by integrating out the other fixed effects except the
## variance term.
obj <- MakeADFun(data=list(x=x,y=y), random=c("intercept", "slope"),
                 parameters=list(intercept=1, slope=0.5, logsd=-2),
                 DLL='linmod')
test <- nlminb(start=-2, obj$fn)
exp(test$par)
