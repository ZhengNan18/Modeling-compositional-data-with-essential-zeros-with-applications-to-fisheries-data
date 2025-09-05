#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{    
    DATA_MATRIX(theta);
    
    PARAMETER(log_sigmaP1);
    PARAMETER(logit_ar1);
    
    PARAMETER(log_sigmaP2);
    PARAMETER(logit_ar2);
    
    PARAMETER_VECTOR(mu);
    
    vector<Type> temp1 = theta.row(0);
    vector<Type> temp2 = theta.col(0);
    
    int n_year = temp2.size();
    int n_age = temp1.size()+1;
    
    
    // matrix<Type> epsilon(n_year,(n_age-1));
    vector<Type> log_theta((n_age-1));
    
    vector<Type> del(n_age-1);
    
    Type one = 1.0;
    //Type zero = 0.0;
    
    using namespace density;
    Type nll = 0.0;
    
    Type sigmaP1 = exp(log_sigmaP1);
    Type ar1 = 2*exp(logit_ar1)/(one + exp(logit_ar1))-1;  // -1<ar<1
    
    Type sigmaP2 = exp(log_sigmaP2);
    Type ar2 = 2*exp(logit_ar2)/(one + exp(logit_ar2))-1; // -1<ar2<1
    
 //   for (int j=0; j<n_year; ++j){
 //       log_theta = theta.row(j);
 //       log_theta = log(log_theta);
 //       nll -= dnorm(log_theta, mu, sigmaP2, true).sum()-sum(log_theta);   //theta ~ log normal
 //   }
    
    del = mu;
    nll += SCALE(AR1(ar1),sigmaP1)(del);
    
    
    for(int i = 0;i < n_year; ++i){
      //vector<Type> del = epsilon.row(i);
      vector<Type> log_theta_temp = theta.row(i);
        log_theta_temp = log(log_theta_temp);
        del = log_theta_temp - mu;
      nll += SCALE(AR1(ar2),sigmaP2)(del);
    }
    
    REPORT(ar1);
    REPORT(ar2);
    REPORT(sigmaP1);
    REPORT(sigmaP2);
    REPORT(mu);
    
    ADREPORT(mu);
  
  return nll;
}

