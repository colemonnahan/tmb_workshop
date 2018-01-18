// TMB implementation of Godwits example

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(pairs);
  DATA_VECTOR(cov1);
  DATA_VECTOR(cov2);
  DATA_VECTOR(cov3);
  DATA_FACTOR(yr); 		// index for year
  // Parameters
  PARAMETER(beta0);		// global intercept
  PARAMETER(beta1);
  PARAMETER(beta2);
  PARAMETER(beta3);
  PARAMETER_VECTOR(tau);
  PARAMETER(logsigma_tau);

  Type sigma_tau = exp(logsigma_tau);
  int n = pairs.size();

  // predictor
  vector<Type> pred(n);
  Type linpred;
  for(int i=0; i<n; i++){
    // Prediction on the linear scale
    linpred = beta0+beta1*cov1(i) + beta2*cov2(i) + beta3*cov3(i) +
      tau(yr(i));
    // log link for Poisson keeps the expected value positive
    pred(i)=exp(linpred);
  }

  Type nll = -dpois(pairs, pred, true).sum();
  nll -= dnorm(tau, 0, sigma_tau, true).sum();
  // Reporting
  REPORT(pred);
  return nll;
}
