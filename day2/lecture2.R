## Lecture 2 code

## Variable types in R
x <- list(a=5L, b=5, c=rnorm(5), d='hello',
          e=factor(c('low', 'high')))
str(x)

## Look at big_orange model
inputs <- readRDS("tmb_models/orange_inputs.RDS")
data <- inputs$data
pars <- inputs$pars
str(data)
str(pars)

## Demonstrate vector behavior in R
pred <- 1:5
pred[6]
pred[6] <- 6
pred
## pred <- 1:5
## > pred[6] [1] NA
## > pred[6] <- 6
## > pred [1] 1 2 3 4 5 6

### Demonstration of methods of uncertainty
## Look at NLL shapes for contrived Poisson
nll <- function(lam, n) sapply(lam, function(i) -sum(dpois(x=rep(4,n), i, TRUE)))
lam <- seq(1, 8, len=1000)
y1 <- nll(lam, 2)
y2 <- nll(lam, 8)
y3 <- nll(lam, 20)
png('likelihoods.png', width=3, height=2, units='in', res=300)
par(mfrow=c(1,1), mar=c(1,2,.5,.5), mgp= c(1, .05, 0), cex.axis=.7, tck=-0.01)
plot(lam, y1-min(y1), type='l', lwd=2, xlab='lambda', ylab='NLL-min(NLL)')
lines(lam, y2-min(y2), col=2, lty=2, lwd=2)
lines(lam, y3-min(y3), col=4, lty=3, lwd=2)
dev.off()
gr <- function(lam, n){
  h <- .00001
  d <- (nll(lam+h, n=n)-nll(lam, n=n))/h
  return(d)
}
## Make plot and table of hessian and SEs
gr1 <- sapply(lam, function(ll) gr(ll, n=2))
gr2 <- sapply(lam, function(ll) gr(ll, n=8))
gr3 <- sapply(lam, function(ll) gr(ll, n=20))
he1 <- optim(4, nll, n=2, hessian=TRUE)$hessian
he2 <- optim(4, nll, n=8, hessian=TRUE)$hessian
he3 <- optim(4, nll, n=20, hessian=TRUE)$hessian
he <- c(he1, he2, he3)
png('derivs.png', width=3, height=2, units='in', res=300)
par(mfrow=c(1,1), mar=c(1,2,.5,.5), mgp= c(1, .05, 0), cex.axis=.7, tck=-0.01)
plot(lam, gr1, xlim=c(3,5), type='l', ylim=c(-5,2), ylab='1st Deriv', xlab=NA)
lines(lam, gr2, col=2, lty=2)
lines(lam, gr3, col=4, lty=3)
plot_derv <- function(g, x=4, eps=.15, col=2){
  y <- 0
  points(x,y,pch=16)
  xx <- c(x-eps, x+eps)
  yy <- c(y-eps*g, y+eps*g)
  lines(xx,yy, col=col, lwd=2.5)
}
plot_derv(he1, col=1)
plot_derv(he2, col=2)
plot_derv(he3, col=4)
cbind(he, 1/sqrt(he))
dev.off()

## 2d NLL shapes
x1 <- x2 <- seq(-2,2, len=25)
x <- expand.grid(x1,x2)
z <- sapply(1:nrow(x), function(i) -sum(dnorm(as.numeric(x[i,]), 0, 1, TRUE)))
z1 <- matrix(z, nrow=25, ncol=25, byrow=TRUE)
z <- sapply(1:nrow(x), function(i) -sum(dnorm(as.numeric(x[i,]), 0, .6, TRUE)))
z2 <- matrix(z, nrow=25, ncol=25, byrow=TRUE)
zlim <- range(c(z1, z2))
par(mfrow=c(2,2))
persp(x1, x2, z1, xlim=c(-2,2), ylim=c(-2,2), zlim=zlim, col='lightblue')
persp(x1, x2, z2, xlim=c(-2,2), ylim=c(-2,2), zlim=zlim, col='lightblue')
contour(x1,x2,z1, xlim=c(-2,2), ylim=c(-2,2), zlim=zlim)
contour(x1,x2,z2, xlim=c(-2,2), ylim=c(-2,2), zlim=zlim)

## A confounded model to demonstrate failure of Hessian inversion
x1 <- rnorm(50); x2 <- x1; y <- 1+x1+rnorm(50)
lm(y~x1+x2)
mod <- lm(y~x1+x2)
summary(mod)

## Working with asymptotic uncertainties
compile("tmb_models/bevholt2.cpp")
dyn.load(dynlib("tmb_models/bevholt2"))
dat <- read.table("tmb_models/bevholt.dat", header=TRUE)
data <- list(SSB=dat$ssb,logR=dat$logR)
parameters <- list(logA=0, logB=0, logsigma=0)
obj <- MakeADFun(data,parameters,DLL="bevholt2")
obj$env$beSilent() # silences console output
opt <- with(obj, nlminb(par, fn, gr))
mles <- opt$par
rep <- sdreport(obj)
summary(rep, 'fixed')
summary(rep, 'report')
fit.cov  <- rep$cov.fixed
ses <- sqrt(diag(fit.cov))
## Delta method, explicitly do y=f(x) is y=exp(logA)
## Var(f(x))=Var(x)*[f'(x)]^2
## Var(f(logA))=Var(exp(logA))=Var(logA)*[exp(logA)]^2
ses[1]*(exp(mles[1]))
## Compare to TMB output (should match)
rep$sd[1]

## Likelihood profile
profA <- tmbprofile(obj, 'logA', parm.range=c(1.5, 2.5))
profB <- tmbprofile(obj, 'logB', parm.range=c(-14, -10))
plot(profA)
confint(profA)

## Bootstrapping demo
Nboot <- 1000
boot.results <- matrix(NA, nrow=Nboot, ncol=3)
N <- length(data$SSB)
for(i in 1:Nboot){
  ind <- sample(1:N, size=N, replace=TRUE)
  data2 <- list(SSB=data$SSB[ind], logR=data$logR[ind])
  obj <- MakeADFun(data2,parameters,DLL="bevholt2")
  obj$env$beSilent() # silences console output
  opt.temp <- with(obj, nlminb(par, fn, gr))
  boot.results[i,] <- opt.temp$par
}

## Plot them together
dev.off()
(CI.A <- as.vector(confint(profA)))
(CI.B <- as.vector(confint(profB)))
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(1.5,.5,.01))
hist(boot.results[,1], main='logA', freq=FALSE, xlab=NA)
box()
xx <- seq(1,3, len=1000)
lines(xx, dnorm(xx, mean=rep$par.fixed[1], sd=ses[1]), lwd=2)
lines(c(CI.A), y=c(1.5,1.5), lwd=2, lty=2)
hist(boot.results[,2], main='logB', freq=FALSE, xlab=NA)
box()
xx <- seq(-14, -11, len=1000)
lines(xx, dnorm(xx, mean=rep$par.fixed[2], sd=ses[2]), lwd=2)
lines(c(CI.B), y=c(.5,.5), lwd=2, lty=2)
quantile(boot.results[,1], probs=c(0.025, 0.975))



