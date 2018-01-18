## Peregrin falcon example

## Read in and prepare the data for modeling.
peregrin <- read.table("tmb_models/falcons.txt", header = TRUE)
yr <- factor(peregrin$Year)         # Create a factor year
mny <- mean(peregrin$Year)
sdy <- sd(peregrin$Year)
cov1 <- (peregrin$Year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
pairs <- peregrin$Pairs

## Fit with lme4
library(lme4)
glmm.fit <- glmer(pairs ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson)
glmm.fit

library(glmmTMB)
glmm.tmbfit <- glmmTMB(pairs ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson)
glmm.tmbfit

## Recreate in TMB
library(TMB)
dat <- list(pairs=pairs, cov1=cov1, cov2=cov2, cov3=cov3, yr=yr)
pars <- list(beta0=4, beta1=0, beta2=0, beta3=0, tau=rep(0,40), logsigma_tau=0)
compile( "tmb_models/falcons.cpp" )
dyn.load( dynlib("tmb_models/falcons") )
obj = MakeADFun(data=dat, parameters=pars, random='tau', DLL="falcons")
obj$fn()
obj$env$beSilent()
opt <- with(obj, nlminb(par, fn, gr))
opt$par
## Predictions from TMB
pred <- obj$report()$pred
## Predictions from lme4.
pred.lme4 <- predict(glmm.fit, type='response')
head(pred-pred.lme4)
# Plot fit
plot(peregrin$Year, pairs, type='b', ylab = "Population size",
     xlab = "Year")
lines(peregrin$Year, pred, lwd = 3, col = "red")

## Compare AICs
AIC(glmm.fit)
AIC(glmm.tmbfit)
2*opt$objective +2*length(opt$par)

## Exercise: Refit without cov2 or cov3 and plot to data
glmm.fit2 <- glmer(pairs ~ (1 | yr) + cov1, family = poisson)
map <- list(beta2=factor(NA), beta3=factor(NA))
obj2 = MakeADFun(data=dat, parameters=pars, random='tau', DLL="falcons", map=map)
obj2$env$beSilent()
opt2 <- with(obj2, nlminb(par, fn, gr))
## They still match glmer
opt2$par
summary(glmm.fit2)

## Calculate AIC
AIC1 <- 2*opt$objective +2*length(opt$par)
AIC2 <- 2*opt2$objective +2*length(opt2$par)
c(AIC1, AIC2)

## Predictions from TMB
pred <- obj2$report()$pred
lines(peregrin$Year, pred, lwd = 3, col = "green")
tt <- as.vector(exp(fixef(glmm.fit2)[1] + fixef(glmm.fit2)[2]*cov1 + unlist(ranef(glmm.fit2))))
lines(peregrin$Year, tt, col='blue')

## Use REML with TMB. We do this by delcaring *everything* but the variance
## as random. That is, all of the fixed effects.
map <- list(beta2=factor(NA), beta3=factor(NA))
obj3 = MakeADFun(data=dat, parameters=pars,
                 random=c('tau','beta1', 'beta0'),
                 DLL="falcons", map=map)
obj3$env$beSilent()
obj3$gr()
opt3 <- with(obj3, nlminb(par, fn, gr))
## Compare estimates of variances. Notice that the REML version is
## bigger. This is because ML estimates of variances are underestimated
## (too small). I don't think REML is possible with glmer or other GLMM
## packages except for TMB.
exp(opt3$par)                           # with REML
exp(opt2$par[3])                        # with ML

