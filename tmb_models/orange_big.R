## Slightly modified version of the scaled up version of the Orange Tree
## example (5,000 latent random variables).

## Demonstration of TMB workflow

## First step: compile model
library(TMB)
compile("orange_big.cpp")
dyn.load(dynlib("orange_big"))

## Second step: build the TMB "object" by passing parameters and data
inputs <- readRDS("orange_inputs.RDS")
obj <- MakeADFun(data=inputs$data, parameters=inputs$pars,
                 random=c("u"), DLL="orange_big")

## Fit model (find minimum negative log likelihood)
opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr)

## Get standard errors of parameters
rep <- sdreport(obj)
rep
