#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{    
  DATA_MATRIX(Po);
  DATA_IVECTOR(count_not0);
  DATA_MATRIX(Po_ratio);
  DATA_IVECTOR(last_not0);   
  DATA_VECTOR(Xop);   
  
//  PARAMETER(log_sigmaP);
  PARAMETER(log_lambda0); 
  PARAMETER(logit_ar0);    
  PARAMETER(log_sd0);
  PARAMETER(logit_ar);    
  PARAMETER(log_sd);
  PARAMETER(log_phi);
  PARAMETER_VECTOR(theta);   
  PARAMETER_VECTOR(dev_p0);   
  
  int n_year = Xop.size(); 
  int n_age = theta.size()+1;  
  Type one = 1.0;
  Type zero = 0.0;
         
//  Type sigmaP = exp(log_sigmaP);
  Type lambda0 = exp(log_lambda0); 
  Type ar0 = 2/(one + exp(-logit_ar0))-1;   
  Type sd0 = exp(log_sd0); 

  Type ar = 2/(one + exp(-logit_ar))-1;   
  Type sd = exp(log_sd); 
  Type phi = exp(log_phi); 

  vector<Type> etheta(n_age);
  
  for(int i = 0;i < n_age-1;++i){etheta(i)=exp(theta(i));}
  etheta(n_age-1) = 1.0;
  vector<Type> p = etheta/sum(etheta);  
  vector<Type> alpha = phi*p;  

  vector<Type> log_p0=log_lambda0+dev_p0;
  vector<Type> p0 = exp(log_p0);
  vector<Type> log_p1 = -log(1+p0);  
  log_p0 = log_p0 + log_p1;
  p0 = exp(log_p0); 
      
  using namespace density;
  Type nll = 0.0;  
    
// Catch age composition NB nll;  
  Type theta_last;
  Type theta_now;
  Type alpha0_now;

  for(int i = 0;i < n_year;++i){
    alpha0_now = 0;
    for(int j = 0;j < n_age;++j){
      if(Po(i,j)==0){
        nll -= log_p0(j); 
      }else{
        nll -= log_p1(j);
        alpha0_now += alpha(j);
        nll += lgamma(alpha(j));
        nll -= (alpha(j)-1)*log(Po(i,j));
      }
    }
    nll -= lgamma(alpha0_now);
  }   

//catch age comp over-dispersion effect nll;

  nll += SCALE(AR1(ar),sd)(theta);
  nll += SCALE(AR1(ar0),sd0)(dev_p0);
  
// predicted residuals;  

  vector<Type> p_true = (1-p0)*p; 
  p_true = p_true/sum(p_true); 
  
  REPORT(p); 
  REPORT(etheta);
  REPORT(ar);
  REPORT(lambda0);
  REPORT(p_true);
  REPORT(sd); 
  REPORT(phi); 
  
  ADREPORT(p);   
  ADREPORT(p_true);  
  
  return nll;
}

