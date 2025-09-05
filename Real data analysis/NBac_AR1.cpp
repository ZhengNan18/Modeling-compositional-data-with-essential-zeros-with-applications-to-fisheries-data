#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{

// input data;  
  DATA_MATRIX(Po);
  DATA_VECTOR(Xop);
  DATA_SCALAR(upbound);   
  
  PARAMETER(log_sigmaP);
  PARAMETER(logit_ar);
  PARAMETER(log_sd_theta);
  PARAMETER(logit_ar_theta);
  PARAMETER_VECTOR(theta);   
  PARAMETER_MATRIX(log_gam); 
//  PARAMETER_MATRIX(log_gamp);  
  // PARAMETER_VECTOR(log_gamt);
  
  int n_year = Xop.size(); 
  int n_age = theta.size()+1;  
  Type one = 1.0;
  int i,j;  
         
  Type sigmaP = exp(log_sigmaP); 
  Type ar = upbound*exp(logit_ar)/(one + exp(logit_ar)); 

  Type sd_theta = exp(log_sd_theta);
  Type ar_theta = 2/(1+exp(-logit_ar_theta))-1;
    
  vector<Type> etheta(n_age);
  vector<Type> gi(n_age);
  vector<Type> gi2(n_age);
    
  vector<Type> pi(n_age); 
  vector<Type> xi(n_age);
  vector<Type> mui(n_age);
  matrix<Type> rpm(n_year,n_age);  
  matrix<Type> resid(n_year,n_age);  
  matrix<Type> residp(n_year,n_age);
  matrix<Type> std_resid(n_year,n_age); 
  matrix<Type> std_residp(n_year,n_age);
  
  using namespace density;
  Type nll = 0.0;  
  
  for(i = 0;i < n_age-1;++i){etheta(i)=exp(theta(i));}
  etheta(n_age-1) = 1.0;
  vector<Type> p = etheta/sum(etheta);
    
// Catch age composition NBac nll;  
  
  for(i = 0;i < n_year;++i){
    xi = Xop(i)*Po.row(i);
    gi = exp(vector<Type>(log_gam.row(i))); 
    pi = p*gi;  
    pi = pi/sum(pi);
    nll -= dmultinom(xi,pi,true);
    rpm.row(i) = pi; 
    resid.row(i) = vector<Type>(Po.row(i)) - pi; 
    std_resid.row(i) = vector<Type>(resid.row(i))/sqrt(pi*(one-pi)/Xop(i));
  }

//catch age comp over-dispersion effect nll;
  for(i = 0;i < n_year;++i){
    vector<Type> del = log_gam.row(i);
    nll += SCALE(AR1(ar),sigmaP)(del); 
    
    //del = log_gamp.row(i);
    //nll += SCALE(AR1(ar),sigmaP)(del);
  }   

  nll += SCALE(AR1(ar_theta),sd_theta)(theta);
     
      
       // simulation code ========================
         SIMULATE {
           array<Type> gi_array(n_age,10000);
           vector<Type> p_true(n_age);
           p_true = 0;
           for(i = 0;i < 10000;++i){
              SCALE(AR1(ar),sigmaP).simulate(gi);
              gi = exp(gi);
              gi_array.col(i) = gi;
              pi = p*gi;
              p_true += pi/sum(pi);
           }
           p_true = p_true/10000;
           REPORT(p_true);
           REPORT(gi_array);
         }

       //==========================================
       
  vector<Type> logit_p = log(p/(1-p));
  vector<Type> p_unbiased = p*exp(sigmaP*sigmaP/2);
  p_unbiased = p_unbiased/sum(p_unbiased);  

  vector<Type> logit_p_true = log(p_unbiased/(1-p_unbiased));

 
                  
   
  REPORT(p); 
  REPORT(rpm);
  REPORT(etheta);
  REPORT(sigmaP);
  REPORT(ar);    
  REPORT(log_gam); 
  REPORT(logit_p); 
  REPORT(logit_p_true); 
  REPORT(p_unbiased); 

  ADREPORT(p);    
  ADREPORT(logit_p);    
  ADREPORT(logit_p_true);    
    
  return nll;
}

