// Simple linear model with slope and interecept
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x)
  DATA_VECTOR(y);
  PARAMETER(intercept);
  PARAMETER(slope);
  PARAMETER(logsd);
  Type sd=exp(logsd);
  Type nll;
  vector<Type> mu(y.size());
  mu=intercept+slope*x;
  nll=-sum(dnorm(y,mu,sd,true));
  return nll;
}
