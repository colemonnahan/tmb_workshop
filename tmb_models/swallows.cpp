#include <TMB.hpp>

// Hand-coded Cauchy distribution
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type inv_logit(Type x){
  Type y= 1/(1+exp(-x));
  return(y);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(I);
  DATA_INTEGER(K);
  DATA_INTEGER(nfam);
  DATA_MATRIX(CH);
  DATA_VECTOR(carez);
  DATA_IVECTOR(year);
  DATA_VECTOR(agec);
  DATA_IVECTOR(family);
  DATA_IVECTOR(last);

  // fixed effects -- bounds added in R
  PARAMETER(sigmayearphi);
  PARAMETER(sigmaphi);
  PARAMETER(sigmap);
  PARAMETER_VECTOR(a);
  PARAMETER(a1);
  PARAMETER_VECTOR(b0);
  PARAMETER_VECTOR(b1);
  // non-centered random effects
  PARAMETER_VECTOR(fameffphi_raw);
  PARAMETER_VECTOR(fameffp_raw);
  PARAMETER_VECTOR(yeareffphi_raw);

  Type nlp=0.0; // negative log prior
  Type nll=0.0; // negative log likelihood
  matrix<Type> p(I,K);
  matrix<Type> phi(I,K-1);
  matrix<Type> chi(I,K+1);

  // To bound below by zero
  Type sigmayearphi2=exp(sigmayearphi);
  Type sigmaphi2=exp(sigmaphi);
  Type sigmap2=exp(sigmap);
  // Jacobian adjustment for variances
  nll -= sigmaphi + sigmayearphi + sigmap;

  p.setZero();
  phi.setZero();
  chi.setZero();

  int k;
  Type x;
  // TMB indexes from 0 not 1, so need to be careful to adjust that
  // below. I've added (-1) where needed.
  for(int i=0; i<I; i++){ // loop over each individual
    // calculate phi as a function of fixed and random effects
    for(int t=0; t<(K-1); t++) {
      x=a(t)+ a1*carez(i)+
      	sigmayearphi2*yeareffphi_raw(year(i)-1)+
      	sigmaphi2*fameffphi_raw(family(i)-1);
      phi(i,t) = inv_logit(x);
    }
    // calculate p as a function of fixed and random effects
    p(i,1-1) = 1;  // first occasion is marking occasion
    for(int t=1; t<K; t++){
      x=b0(year(i)-1)+ b1(year(i)-1)*agec(t)+
      	sigmap2*fameffp_raw(family(i)-1);
      p(i,t) = inv_logit(x);
    }
    // probabilitiy of never being seen after last observation. ind here is
    // a reverse index so this loop goes from K:2, recursively calculating
    // backward.
    chi(i,K+1-1) = 1.0;
    k = K;
    while (k > 1) {
      chi(i,k-1) = (1 - phi(i,k-1-1)) + phi(i,k-1-1) * (1 - p(i,k-1)) * chi(i,k+1-1);
      k = k - 1;
    }
    chi(i,1-1) = (1 - p(i,1-1)) * chi(i,2-1);
  }

  // priors
  nlp-= dnorm(b0, Type(0.0), Type(5), true).sum();
  nlp-= dnorm(b1, Type(0.0), Type(5), true).sum();
  nlp-= dnorm(a, Type(0.0), Type(1.5), true).sum();
  nlp-= dnorm(a1, Type(0.0), Type(5), true);
  nlp-= dcauchy(sigmaphi2, Type(0), Type(1.0), true);
  nlp-= dnorm(sigmayearphi2, Type(0), Type(0.5), true);
  nlp-= dcauchy(sigmap2, Type(0), Type(1.0), true);

  // random effects; non-centered
  nll-=dnorm(fameffphi_raw, Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(fameffp_raw,Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(yeareffphi_raw, Type(0.0), Type(1.0), true).sum();

  // // likelihood
  for(int i=0; i<I; i++){ // loop over each individual
    // probability of survival, known alive since k<last
    for (int t=1; t<last(i); t++) {
    	nll-= log(phi(i,t-1));
    }
    // // probability of observation given known alive
    for(int t=0; t< last(i); t++){
      if(CH(i,t)==1){
    	nll-= log(p(i,t));
      } else {
    	nll-= log(1-p(i,t));
      }
    }
    // probability of no observations after time period last
     nll-= log(chi(i,last(i)+1-1));
  }
  Type nld=nll+nlp; // negative log density
  REPORT(p);
  REPORT(chi);
  REPORT(phi);
  return(nld);
}
