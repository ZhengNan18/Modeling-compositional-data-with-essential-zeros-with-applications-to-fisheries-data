#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data;
  DATA_MATRIX(Po);         // matrix of proportions
  DATA_VECTOR(Xop);        // total counts per year
  
  // Parameters;
  PARAMETER(logit_ar);          
  PARAMETER(log_sigmaP);            
  PARAMETER_VECTOR(theta);     // n_year x (n_age - 1)
  
  // Dimensions;
  int n_year = Xop.size(); 
  int n_age = theta.size() + 1;
  Type one = 1.0;
  
  // Transformed parameters;
  Type ar = 2*exp(logit_ar)/(one + exp(logit_ar))-1; 
  Type sigmaP = exp(log_sigmaP);
  
  
  vector<Type> etheta = exp(theta);        // element-wise exponentiation
  Type denom = 1.0 + etheta.sum();         // denominator for additive logistic
  
  vector<Type> p(n_age);                   // output probability vector
  
  for (int i = 0; i < n_age - 1; ++i) {
    p(i) = etheta(i) / denom;
  }
  p(n_age - 1) = 1.0 / denom;              // fixed last category
  
  
  // Negative log-likelihood
  Type nll = 0.0; 
  
  for(int i = 0; i < n_year; ++i){
    vector<Type> xi = Xop(i) * Po.row(i);      // multinomial counts
    nll -= dmultinom(xi, p, true);            // log-likelihood
  }
  
  // AR(1) penalty on theta across years (by age class)
  using namespace density;
  nll += SCALE(AR1(ar), sigmaP)(theta);
  
  // Optional: If theta row penalty (across ages) is also desired (be careful of double-penalizing)
  // for(int i = 0; i < n_year; ++i){
  //   vector<Type> del = theta.row(i);
  //   nll += SCALE(AR1(ar), sd)(del); 
  // }
  
  // Outputs
  REPORT(theta);
  REPORT(ar);
  REPORT(sigmaP);
  REPORT(p);
  ADREPORT(p)
  
  return nll;
}
