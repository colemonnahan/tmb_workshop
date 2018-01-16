// Beverton-Hold stock recruit function. Modified from A. Nielson's
// lectures: http://www.nielsensweb.org/bergen/

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(SSB);
  DATA_VECTOR(logR);

  PARAMETER(logA);
  PARAMETER(logB);
  PARAMETER(logsigma);

  // We use log transforms to keep parameters positive
  Type B=exp(logB);
  Type sigma=exp(logsigma);

  // Vectorized calculation of prediction vector
  vector<Type> pred= logA+log(SSB)-log(Type(1)+B*SSB);
  // negative log-likelihood
  // Note: using vector calculations so need to sum them
  Type nll=-dnorm(logR,pred, sigma, true).sum();

  // Reporting
  Type A=exp(logA);
  ADREPORT(A);
  ADREPORT(B);

  return nll;
}
