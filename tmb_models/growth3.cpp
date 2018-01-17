#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nobs);	   // number of total observations (rows)
  DATA_VECTOR(loglengths); // observed lengths in log space
  DATA_FACTOR(fish);	   // factor relating observation to fish
  DATA_FACTOR(ages);	   // independent variable
  DATA_VECTOR(ages_pred);  // ages to predict Length at

  // fixed effects
  PARAMETER(logsigma_obs);	// observation variance
  
  // random effect vectors
  PARAMETER_VECTOR(logLinf);	// L infinity in log space
  PARAMETER_VECTOR(logk);	// growth k in log space

  // hyperparameters
  PARAMETER(logsigma_logLinf);	// variance of logLinf (in log space)
  PARAMETER(mean_logLinf);	// mean of logLinf (in log space)
  PARAMETER(logsigma_logk);	// variance of logk (in log space)
  PARAMETER(mean_logk);		// mean of logk (in log space)

  // Exponentiate variances
  Type sigma_logLinf=exp(logsigma_logLinf);
  Type sigma_logk=exp(logsigma_logk);
  Type sigma_obs = exp(logsigma_obs);

  // Calculate predicted lengths in log space
  vector<Type> ypred(Nobs);
  Type k;
  for(int i=0; i<Nobs; i++){
    k = exp(logk(fish(i)-1));
    // t0 assumed known at 5
    ypred(i) = logLinf(fish(i)-1) + log(1-exp(-k*(ages(i)-5)));
  }

  Type nll=0.0; // negative log likelihood initialized at 0
  // likelihood of data
  nll-=dnorm(ypred, loglengths, sigma_obs, true).sum();
  // probabilities of the random effects
  nll-=dnorm(logk, mean_logk, sigma_logk, true).sum();
  nll-=dnorm(logLinf, mean_logLinf, sigma_logLinf, true).sum();

  ADREPORT(logLinf);
  ADREPORT(logk);
  return(nll);
}
