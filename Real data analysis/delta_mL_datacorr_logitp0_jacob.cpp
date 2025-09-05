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
  PARAMETER_VECTOR(logit_ar);    
  PARAMETER(log_sd);
  PARAMETER(logit_ar_data);    
  PARAMETER(log_sd_data);

  PARAMETER(p_mean);

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

  vector<Type> ar = 2/(one + exp(-logit_ar))-1;   
  Type sd = exp(log_sd); 
  Type ar_data = 2/(one + exp(-logit_ar_data))-1;   
  Type sd_data = exp(log_sd_data); 

  matrix<Type> corr_mat(n_age-1,n_age-1); 
  corr_mat.setZero();

  Type corr_now;
  for(int i = 0;i < n_age-1;++i){
    corr_mat(i,i) = sd_data;
    corr_now = sd_data;
    for(int j = i+1;j < n_age-1;++j){
      corr_now = corr_now*ar_data;
      corr_mat(i,j) = corr_now;
      corr_mat(j,i) = corr_now;
    }
  }
    
  vector<Type> etheta(n_age);
  vector<Type> gi(n_age);
  vector<Type> pi(n_age); 
  vector<Type> xi(n_age); 
  vector<Type> mui(n_age);  
  matrix<Type> mu(n_year,n_age); 
  matrix<Type> resid(n_year,n_age); 
  matrix<Type> residp(n_year,n_age);
  matrix<Type> std_resid(n_year,n_age); 
  matrix<Type> std_residp(n_year,n_age);   
  
  for(int i = 0;i < n_age-1;++i){etheta(i)=exp(theta(i));}
  etheta(n_age-1) = 1.0;
  vector<Type> p = etheta/sum(etheta);  

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
  
  for(int i = 0;i < n_year;++i){
    vector<Type> del_now(count_not0(i)-1);
    matrix<Type> jacob_now(count_not0(i)-1,n_age-1); 
    jacob_now.setZero();
    int i_not0 = 0;
    if(last_not0(i)<(n_age-1)){
      theta_last = theta(last_not0(i));
    }else{
      theta_last = 0.0;
    }
    for(int j = 0;j < n_age;++j){
      if(Po_ratio(i,j)==0){
        nll -= log_p0(j); 
      }else{
        nll -= log_p1(j);
        nll += log(Po(i,j));
        if(j<last_not0(i)){
          del_now(i_not0) = log(Po_ratio(i,j))-theta(j)+theta_last;
          jacob_now(i_not0,j) = 1.0;
          if(last_not0(i)<(n_age-1)){
            jacob_now(i_not0,last_not0(i)) = -1.0;
          }
          i_not0 += 1;
        }
      }
    }
    matrix<Type> cov_now = jacob_now*corr_mat*jacob_now.transpose(); 
    nll += MVNORM(cov_now)(del_now);
  }   

//catch age comp over-dispersion effect nll;

  nll += SCALE(AR1(ar(0)),sd)(theta-p_mean);
  nll += SCALE(AR1(ar0),sd0)(dev_p0);
  
// predicted residuals;  

  vector<Type> p_true = (1-p0)*p; 
  p_true = p_true/sum(p_true); 
  vector<Type> logit_p_true = log(p_true/(1-p_true));

  REPORT(p); 
  REPORT(etheta);
  REPORT(ar);
  REPORT(lambda0);
  REPORT(p_true);
  REPORT(sd); 
  REPORT(sd_data); 
  REPORT(ar_data); 
  REPORT(logit_p_true);
  REPORT(theta);
  REPORT(p0);
  
  ADREPORT(p);   
  ADREPORT(p_true);  
  ADREPORT(logit_p_true);
  
  return nll;
}

