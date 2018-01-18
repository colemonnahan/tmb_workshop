#include <TMB.hpp>

//// ---------- Custom likelihood functions, used be used in template
//// below. These are not built into TMB like dnorm and dgamma are.
//log-normal likelihood
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse Gaussian
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}

//// ---------- The CPUE model
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(likelihood); 	// Likelihood flag
  DATA_VECTOR(y);		// Observed catches (can be zero)
  DATA_VECTOR(lat);		// latitude
  DATA_VECTOR(lon);		// longitude
  DATA_IVECTOR(out_sample)		// boolean for out-of-sample or not

  // Parameters
  PARAMETER(intercept);		// intercept
  PARAMETER(theta); 		// logit of zero_prob
  PARAMETER(beta_lat);		// slope for longtitude
  PARAMETER(beta_lon);		// slope for latititude
  PARAMETER(logsigma);		// log of observation SD

  // Objective funcction
  Type zero_prob = 1 / (1 + exp(-theta));
  Type sigma = exp(logsigma);
  int n = y.size();
  Type jnll=0;	 // used jnll
  Type pred_jnll=0;		// predicted NLL (out of sample)

  // Linear predictor
  vector<Type> pred(n);
  pred = intercept + beta_lat*lat + beta_lon*lon;

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n; i++){
    Type jnlltemp=0; // NLL for this data point
    // If data is zero
    if(y(i)==0){
      jnlltemp -= log( zero_prob );
    } else {
      // Positive data contribution
      jnlltemp -= log( 1-zero_prob );
      // And the contribution of the likelihood
      if(likelihood==1)
	jnlltemp -= dinvgauss(y(i), exp(pred(i)), sigma, true);
      else if(likelihood==2)
	jnlltemp -= dlognorm(y(i), pred(i), sigma, true);
      else {
	std::cout << "Invalid likelihood specified" << std::endl;
	return 0;
      }
    }
    if( out_sample(i)==0 ) jnll += jnlltemp;
    if( out_sample(i)==1 ) pred_jnll += jnlltemp;
  }

  // Reporting
  REPORT(zero_prob);
  REPORT(pred);
  REPORT(sigma);
  REPORT(pred_jnll);
  return jnll;
}
