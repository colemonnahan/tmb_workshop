// Simulated von bertalanffy growth with single k, Linf
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nfish);	   // number of fish
  DATA_INTEGER(Nobs);	   // number of total observations (rows)
  DATA_VECTOR(loglengths); // observed lengths in log space
  DATA_FACTOR(fish);	   // factor relating observation to fish
  DATA_FACTOR(ages);	   // independent variable
  DATA_VECTOR(ages_pred);  // ages to predict Length at

  // fixed effects
  PARAMETER(logsigma);		// observation variance
  PARAMETER(logLinf);		// L infinity in log space
  PARAMETER(logk);		// growth k in log space

  // Calculate predicted lengths in log space
  vector<Type> ypred(Nobs);
  Type k = exp(logk);
  for(int i=0; i<Nobs; i++){
    // t0 assumed known at 5
    ypred(i) = logLinf + log(1-exp(-k*(ages(i)-5)));
  }

  // likelihood of data
  Type sigma = exp(logsigma);
  Type nll= -dnorm(ypred, loglengths, sigma, true).sum();

  // Calculate lengths (not loglengths!) for inpute ages_pred.
  int Npred = ages_pred.size();
  vector<Type> lengths_pred(Npred);
  for(int i=0; i<Npred; i++){
    lengths_pred(i)=exp(logLinf+log(1-exp(-k*(ages_pred(i)-5))));
  }
  ADREPORT(lengths_pred);
  return(nll);
}
