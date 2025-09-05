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
  DATA_MATRIX(Zo);
  
//  PARAMETER(log_sigmaP);
  PARAMETER(log_lambda0); 
  PARAMETER(logit_ar0);    
  PARAMETER(log_sd0);
  PARAMETER(logit_ar);    
  PARAMETER(log_sd);
//  PARAMETER(log_phia);
  PARAMETER(logit_ar_b);    
  PARAMETER(log_sd_b);
//  PARAMETER(log_phib);
  PARAMETER_VECTOR(lai);   
  PARAMETER_VECTOR(lbi);   
  PARAMETER_VECTOR(dev_p0);   
  
  int n_year = Xop.size(); 
  int n_age = lai.size()+1;  
  Type one = 1.0;
  Type zero = 0.0;
         
  Type lambda0 = exp(log_lambda0); 
  Type ar0 = 2/(one + exp(-logit_ar0))-1;   
  Type sd0 = exp(log_sd0); 

  Type ar = 2/(one + exp(-logit_ar))-1;   
  Type sd = exp(log_sd); 
  //Type phia = exp(log_phia);
  Type ar_b = 2/(one + exp(-logit_ar_b))-1;   
  Type sd_b = exp(log_sd_b); 
  //Type phib = exp(log_phib);

  vector<Type> ai = exp(lai);
  vector<Type> bi = exp(lbi);
  vector<Type> nui(n_age-1); 
  for(int i = 0;i < n_age-2;++i){
     nui(i) = bi(i)-ai(i+1)-bi(i+1);
  }
  nui(n_age-2) = bi(n_age-2)-1;

  vector<Type> log_p0=log_lambda0+dev_p0;
  vector<Type> p0 = exp(log_p0);
  vector<Type> log_p1 = -log(1+p0);  
  log_p0 = log_p0 + log_p1;
  p0 = exp(log_p0); 

  vector<Type> p(n_age); 

  p(0) = (1-p0(0))*ai(0)/( ai(0) + bi(0) ); 
  Type prod_now = 1; 
  for(int i = 1;i < n_age-1;++i){
     prod_now = prod_now*( 1 - (1-p0(i-1))*ai(i-1)/( ai(i-1) + bi(i-1) ) );
     p(i) = prod_now*(1-p0(i))*ai(i)/( ai(i) + bi(i) );
  }
  p(n_age-1) = prod_now*( 1 - (1-p0(n_age-2))*ai(n_age-2)/( ai(n_age-2) + bi(n_age-2) ) );

  using namespace density;
  Type nll = 0.0;  
    
// Catch age composition NB nll;  
  Type theta_last;
  Type theta_now;
  Type alpha0_now;
  int j_now = 0;
  for(int i = 0;i < n_year;++i){
    for(int j = 0;j < n_age;++j){
      if(Po(i,j)==0){
        nll -= log_p0(j); 
      }else{
        nll -= log_p1(j);
        if(j<last_not0(i)){
           nll -= dbeta(Zo(i,j),ai(j),bi(j),true);            
        }
      }
    }
  }   
  
  Type nll_1 = nll;
//catch age comp over-dispersion effect nll;

  nll += SCALE(AR1(ar),sd)(lai);
  nll += SCALE(AR1(ar_b),sd_b)(lbi);
  nll += SCALE(AR1(ar0),sd0)(dev_p0);
  
// predicted residuals;  

  vector<Type> p_true = p; 
  p_true = p_true/sum(p_true); 
  
  REPORT(p); 
  REPORT(ar);
  REPORT(lambda0);
  REPORT(p_true);
  REPORT(sd); 
  REPORT(log_p0); 
  REPORT(log_p1); 
  REPORT(nui); 
  REPORT(nll_1); 
  
  ADREPORT(p);   
  ADREPORT(p_true);  
  
  return nll;
}

