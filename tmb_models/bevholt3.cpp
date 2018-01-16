// Beverton-Hold stock recruit function. Modified from A. Nielson's
// lectures: http://www.nielsensweb.org/bergen/

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(SSB);
  DATA_VECTOR(logR);
  DATA_VECTOR(SSB_pred);

  PARAMETER(logA);
  PARAMETER(logB);
  PARAMETER(logsigma);
    // We use log transforms to keep parameters positive
  Type B=exp(logB);
  Type sigma=exp(logsigma);
  
 // !!!! USE A LOOP HERE FOR GROWTH MODEL !!!!!!
  // Vectorized calculation of prediction vector
  vector<Type> pred= logA+log(SSB)-log(Type(1)+B*SSB);

  // negative log-likelihood
  // Note: using vector calculations so need to sum them
  Type nll=-dnorm(logR,pred,sigma,true).sum();

  // Predict logR at SSB_pred
  vector<Type> logR_pred=
    logA+log(SSB_pred)-log(Type(1)+B*SSB_pred);
  // Get SE for this
  ADREPORT(logR_pred);
  return nll;
}
