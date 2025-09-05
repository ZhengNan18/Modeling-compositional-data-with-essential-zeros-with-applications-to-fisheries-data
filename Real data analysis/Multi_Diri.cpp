#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data;
  DATA_MATRIX(Po);         
  DATA_VECTOR(Xop);       
  
  // Parameters;
  PARAMETER(log_phi);           
  PARAMETER(logit_ar);        
  PARAMETER(log_sigmaP);            
  PARAMETER_VECTOR(theta);      
  
  // Dimensions;
  int n_year = Xop.size(); 
  int n_age = theta.size() + 1;
  
  // Transformed parameters
  Type one = 1.0;
  Type phi = exp(log_phi);
  Type ar = 2*exp(logit_ar)/(one + exp(logit_ar))-1;  
  Type sigmaP = exp(log_sigmaP);
  
  // Unnormalized age composition: etheta = exp(theta), last element fixed to 1;
  
  vector<Type> etheta(n_age);
  
  for(int i = 0; i < n_age - 1; ++i) {
    etheta(i) = exp(theta(i));
  }
  etheta(n_age - 1) = 1.0;
  
  // Normalize to get proportions;
  
  vector<Type> p = etheta / etheta.sum();
  vector<Type> alpha = phi * p;  
  
  // Dirichlet-Multinomial Likelihood
  Type nll = 0.0;
  for(int i = 0; i < n_year; ++i) {
    vector<Type> xi = Xop(i) * Po.row(i);  // Convert proportions to counts
    Type N = xi.sum();                     // Total count
    Type alpha_sum = alpha.sum();         // α₀
    
    nll -= lgamma(alpha_sum);
    nll -= lgamma(N + 1.0);
    nll += lgamma(N + alpha_sum);
    
    for(int j = 0; j < n_age; ++j) {
      nll -= lgamma(xi(j) + alpha(j));
      nll += lgamma(alpha(j));
      nll += lgamma(xi(j) + 1.0);
    }
  }  
 
  // Catch age composition over-dispersion effect;
  using namespace density;
  nll += SCALE(AR1(ar), sigmaP)(theta);
  
  vector<Type> logit_p = log(p/(1-p));

  // Outputs for diagnostics
  REPORT(p);
  REPORT(etheta);
  REPORT(ar);
  REPORT(sigmaP);
  REPORT(phi);
  REPORT(logit_p);
  REPORT(alpha);

  ADREPORT(logit_p);
  ADREPORT(p);   
  ADREPORT(alpha);
  
  return nll;
}
