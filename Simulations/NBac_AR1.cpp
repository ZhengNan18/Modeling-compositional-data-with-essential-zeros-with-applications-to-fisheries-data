#include <TMB.hpp> 
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{

// input data;  
  DATA_MATRIX(Po);
  DATA_VECTOR(Xop);    
  
  PARAMETER(log_sigmaP);
  PARAMETER(logit_ar);
  PARAMETER(log_sd_theta);
  PARAMETER(logit_ar_theta);
  PARAMETER_VECTOR(theta);   
  PARAMETER_MATRIX(log_gam); 
  PARAMETER_MATRIX(log_gamp);  
  // PARAMETER_VECTOR(log_gamt);
  
  int n_year = Xop.size(); 
  int n_age = theta.size()+1;  
  Type one = 1.0;
  int i,j;  
         
  Type sigmaP = exp(log_sigmaP); 
  Type ar = 2*exp(logit_ar)/(one + exp(logit_ar))-1; 

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
    
    del = log_gamp.row(i);
    nll += SCALE(AR1(ar),sigmaP)(del);
  }   

  nll += SCALE(AR1(ar_theta),sd_theta)(theta);
     
  for(i = 0;i < n_year;++i){
    xi = Xop(i)*Po.row(i);
    gi = exp(vector<Type>(log_gamp.row(i))); 
    pi = p*gi;  
    pi = pi/sum(pi);
    residp.row(i) = vector<Type>(Po.row(i)) - pi; 
    std_residp.row(i) = vector<Type>(residp.row(i))/sqrt(pi*(one-pi)/Xop(i));
  }
  
  Type RMS_std_residp = 0.0; 
  for(i = 0;i < n_year;++i){
    for(j = 0;j < n_age;++j){
      RMS_std_residp += std_residp(i,j)*std_residp(i,j);
  }}
  RMS_std_residp = sqrt(RMS_std_residp/(n_age*n_year));          
    
      
       // simulation code ========================
         SIMULATE {
           vector<Type> p_true(n_age);
           p_true = 0;
           for(i = 0;i < 10000;++i){
              SCALE(AR1(ar),sigmaP).simulate(gi);
              gi = exp(gi);
              pi = p*gi;
              p_true += pi/sum(pi);
           }
           p_true = p_true/10000;
           REPORT(p_true);
         }

       //==========================================
       
       
 
                  
   
  REPORT(p); 
  REPORT(rpm);
  REPORT(etheta);
  REPORT(sigmaP);
  REPORT(ar);    
  REPORT(log_gam); 
  REPORT(resid);  
  REPORT(residp); 
  REPORT(std_resid);  
  REPORT(RMS_std_residp);   

  
  ADREPORT(p);    
  ADREPORT(residp);
    
  return nll;
}

