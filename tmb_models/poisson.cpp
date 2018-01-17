// Simulated Poisson hierarchical GLM
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(C);		// observed counts
  DATA_FACTOR(sites);		// index for sites
  PARAMETER_VECTOR(D);		// densities at each site
  PARAMETER(mu);		// sites mean
  PARAMETER(logsigma);		// sites variability
  Type sigma=exp(logsigma);

  // Loop through each observed count and make prediction
  int n=C.size();
  vector<Type> pred(n);
  for(int i=0; i<n; i++){
    // Expected counts are just the random effect for that site
    pred(i)=D(sites(i)-1);
  }

  // Likelihood calculations
  Type nll=0;
  // The data component
  nll -= dpois(C,pred, true).sum();
  // The random effects hyperdistribution
  nll -= dnorm(D, mu, sigma, true).sum();
  return nll;
}
