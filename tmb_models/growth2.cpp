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
  PARAMETER_VECTOR(logLinf);	// L infinity in log space
  PARAMETER_VECTOR(logk);		// growth k in log space

  // Calculate predicted lengths in log space
  vector<Type> ypred(Nobs);
  vector<Type> ktemp(Nobs);
  Type k;
  for(int i=0; i<Nobs; i++){
    k = exp(logk(fish(i)-1));
    ktemp(i)=k;
    // t0 assumed known at 5
    ypred(i) = logLinf(fish(i)-1) + log(1-exp(-k*(ages(i)-5)));
  }

  // likelihood of data
  Type sigma = exp(logsigma);
  Type nll= -dnorm(ypred, loglengths, sigma, true).sum();
  REPORT(ktemp);
  return(nll);
}
