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
  DATA_FACTOR(site);		// Random effect index for observation i
    // SPDE objects from R-INLA
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER(intercept);		// intercept
  PARAMETER(theta); 		// logit of zero_prob
  PARAMETER(beta_lat);		// slope for longtitude
  PARAMETER(beta_lon);		// slope for latititude
  PARAMETER(logsigma);		// log of observation SD
  PARAMETER(logsigma_space);	// spatial variance
  PARAMETER_VECTOR(u);		// spatial random effects
  PARAMETER(logkappa);		// kind of the decorelation

  //// Geospatial derived quantities, given parameters
  Type Range = sqrt(8) / exp( logkappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*logsigma_space) * exp(2*logkappa));
  Eigen::SparseMatrix<Type>
    Q = exp(4*logkappa)*M0 + Type(2.0)*exp(2*logkappa)*M1 + M2;

  Type zero_prob = 1 / (1 + exp(-theta));
  Type sigma = exp(logsigma);
  Type sigma_space = exp(logsigma_space);
  int n = y.size();
  Type jnll=0;	 // used jnll

  // Linear predictor
  vector<Type> pred(n);
  for(int i=0; i<n; i++){
    pred(i) = intercept + beta_lat*lat(i) + beta_lon*lon(i) +
      u(site(i));		// the spatial effect
  }
  // Probability of random effects (hyperdistribution). This replaces neg_log_density (MVN)
  // because it is sparse
  using namespace density;
  // Spatial random effects likelihood
  jnll += SCALE(GMRF(Q), 1/sigma_space)(u); // returns negative already
  
  // Probability of data conditional on fixed effect values
  for(int i=0; i<n; i++){
    // If data is zero
    if(y(i)==0){
      jnll -= log( zero_prob );
    } else {
      // Positive data contribution
      jnll -= log( 1-zero_prob );
      // And the contribution of the likelihood
      if(likelihood==1)
	       jnll -= dinvgauss(y(i), exp(pred(i)), sigma, true);
      else if(likelihood==2)
	       jnll -= dlognorm(y(i), pred(i), sigma, true);
      else {
      	std::cout << "Invalid likelihood specified" << std::endl;
        return 0;
      }
    }
  }

  // Reporting
  REPORT(zero_prob);
  REPORT(pred);
  REPORT(sigma);
  REPORT(u);
  REPORT(Range);
  REPORT(SigmaO);
  return jnll;
}
