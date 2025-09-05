library(TMB) 
library('ggplot2')  
library('gridExtra')
#library('ggcorrplot')

mlogit = function(U){
  n.age = length(U[1,])
  cU = t(apply(U,1,cumsum))[,1:(n.age-1)]
  tU = U[,1:(n.age-1)]
  ret = log(tU) - log(1-cU)
  return(ret)
}         

milogit = function(Ui){
  ret = exp(Ui);
  temp = 1+ret;
  den = t(apply(temp,1,cumprod))
  ret=ret/den
  ret = cbind(ret,1/den[,n.age-1])
  return(ret)
} 

# page 18
#scenario 1
pv=c(45,147,124,44,32,21,20,9,1,15,7)
pv=pv/sum(pv)

#scenario 2
# pv=rep(1,11)
# pv=pv/sum(pv)

#scenario 3
# pv=c(1,4,6,8,10,12,10,8,6,4,1)
# pv=pv/sum(pv)


n.age = length(pv)  
n.year=10  
nsim=500

ar1.phi = 0.8
sigmaP = 1 
mean_tot_X=30



compile("NBac_AR1.cpp")  

dyn.load(dynlib("NBac_AR1"))  
# dyn.load("NBac_AR1") 
#dyn.unload("NBac_AR1"); 

compile("delta_mL_datacorr.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")  

dyn.load(dynlib("delta_mL_datacorr")) 
#dyn.unload("delta_mL_datacorr")

# delta clr
compile("delta_clr_fix_datacorr.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")  

dyn.load(dynlib("delta_clr_fix_datacorr")) 
# dyn.load("delta_clr_fix_datacorr") 
#dyn.unload("delta_clr_fix_datacorr")


#### 
compile("DRHT.cpp")  
dyn.load(dynlib("DRHT")) 

# logitp0
compile("delta_mL_datacorr_logitp0.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_mL_datacorr_logitp0"))
#dyn.unload(dynlib("delta_mL_datacorr_logitp0"))

#Dirichlet
compile("delta_Dirichlet.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_Dirichlet"))
#dyn.unload(dynlib("delta_Dirichlet"))

#Dirichlet
compile("delta_Dirichlet_AR1_0.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_Dirichlet_AR1_0"))
#dyn.unload(dynlib("delta_Dirichlet_AR1_0"))

#generalized Dirichlet

compile("delta_generalized_Dirichlet_AR1_0_zi.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_generalized_Dirichlet_AR1_0_zi"))
#dyn.unload(dynlib("delta_generalized_Dirichlet_AR1_0_zi"))

compile("delta_generalized_Dirichlet.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_generalized_Dirichlet"))
#dyn.unload(dynlib("delta_generalized_Dirichlet"))

compile("Multinomial.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign") 

dyn.load(dynlib("Multinomial"))
#dyn.unload(dynlib("Multinomial"))

compile("Multi_Diri.cpp")  

dyn.load(dynlib("Multi_Diri"))
#dyn.unload(dynlib("Multi_Diri"))

set.seed(100)

# filename1 = "plot1.jpeg"
#filename2 = "plot2.jpeg"

## parameters and bounds for NB;

# containers for simulation results;

sim.nbac = matrix(NA,nrow=nsim,ncol=n.age)
sim.nb = matrix(NA,nrow=nsim,ncol=n.age)

sim.nbac_true = matrix(NA,nrow=nsim,ncol=n.age)
sim.nb_true = matrix(NA,nrow=nsim,ncol=n.age)

sim.nbac_AR1 = matrix(NA,nrow=nsim,ncol=n.age)
sim.nbac_AR1_true = matrix(NA,nrow=nsim,ncol=n.age) 

sim.mn = matrix(NA,nrow=nsim,ncol=n.age)   
sim.mL = matrix(NA,nrow=nsim,ncol=n.age)   

sim.delta_mL_datacorr = matrix(NA,nrow=nsim,ncol=n.age)   
sim.delta_clr_fix_datacorr = matrix(NA,nrow=nsim,ncol=n.age) 

sim.nll = matrix(NA,nrow=nsim,ncol=3)    
sim.resvar = matrix(NA,nrow=nsim,ncol=3) 
sim.respvar = matrix(NA,nrow=nsim,ncol=3)    
var.par = matrix(NA,nrow=nsim,ncol=2)

sim.multinomial = matrix(NA,nrow=nsim,ncol=n.age)

sim.Diri_multinomial = matrix(NA,nrow=nsim,ncol=n.age)

pplot=FALSE

data_save = matrix(NA,nrow=nsim,ncol=n.age)

############################################
p_true_save_clr = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_ml = matrix(NA,nrow=nsim,ncol=n.age)
p_true_save_ml_logitp0 = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_DRHT = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_Dirichlet = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_Dirichlet_0AR1 = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_AR1_0 = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_AR1_0_perturb = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_perturb = matrix(NA,nrow=nsim,ncol=n.age)

set.seed(100)
error_count = 0
iter = 1
while(iter <= nsim){
  
  Xop = rpois(n.year,mean_tot_X) 
  #Xop = rep(mean_tot_X,n.year)  
  
  Xo = matrix(NA,nrow=n.year,ncol=n.age)
  Xr = matrix(NA,nrow=n.year,ncol=n.age)
  
  gam.matrix = matrix(NA,nrow=n.year,ncol=n.age)
  
  for(i in 1:n.year){
    gam.matrix[i,] = exp(sigmaP*arima.sim(list(order=c(1,0,0), ar=ar1.phi), n=n.age)*sqrt(1 - ar1.phi**2))
    #gam.matrix[i,] = exp(sigmaP*rnorm(0, n=n.age))
  }
  
  rpm = matrix(pv,nrow=n.year,ncol=n.age,byrow=T)*gam.matrix
  rpm=rpm/matrix(apply(rpm,1,sum),nrow=n.year,ncol=n.age,byrow=F)   
  
  for(i in 1:n.year){
    Xo[i,] = rmultinom(1, Xop[i], rpm[i,])
    Xr[i,] = Xo[i,]
    ind = Xr[i,]==0
    Xr[i,ind]=0.5     ### page 18
    # Xr[i,ind]=pv[ind]*Xop[i]    ### page 19
  }
  
  Po = Xo/matrix(Xop,nrow=n.year,ncol=n.age,byrow=F)
  
  data_save[iter,] = colMeans(Po)
  ####
  Po_not0 = (Po != 0)
  count_not0 = rowSums(Po_not0)
  last_not0 = rep(NA,dim(Po)[1])
  for (i in 1:dim(Po)[1]) {
    last_not0[i] = max(which(Po_not0[i,]))
  }
  
  if(sum(count_not0<=1)>0){next}
  
  Po_ratio = Po
  for (i in 1:dim(Po)[1]) {
    Po_ratio[i,] = Po[i,]/Po[i,last_not0[i]]
  }
  
  
  
  ###########  do estimation ###################; 
  
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(Po=Po,Xop=Xo.guess)
  
  parameters_NBac_AR1 <- list( 
    log_sigmaP=log(0.2),                  
    logit_ar = 0, 
    log_sd_theta = log(0.2),
    logit_ar_theta = 0,
    theta = rep(0,n.age-1), 
    log_gam=matrix(0,nrow=n.year,ncol=n.age), 
    log_gamp=matrix(0,nrow=n.year,ncol=n.age)
    # log_gamt = rep(0, n.age)
  )
  
  NBac_AR1.obj <- MakeADFun(tmb.data,parameters_NBac_AR1,random=c("theta","log_gam","log_gamp"), 
                            DLL="NBac_AR1",inner.control=list(maxit=5000,trace=F),silent = TRUE)  
  
  NBac_AR1.opt<-nlminb(NBac_AR1.obj$par,NBac_AR1.obj$fn,NBac_AR1.obj$gr,
                       control = list(trace=0,eval.max=2000,iter.max=1000))
  
  NBac_AR1.rep = NBac_AR1.obj$report() 
  NBac_AR1.rep = NBac_AR1.obj$simulate()
  
  sim.nbac_AR1[iter,] = NBac_AR1.rep$p 
  sim.nbac_AR1_true[iter,] = NBac_AR1.rep$p_true 
  ########################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(Po=Po,Xop=Xo.guess)
  
  ## parameters and bounds for NBac;
  
  parameters <- list(
    logit_ar = qlogis(0.5),
    log_sigmaP=log(0.1),
    theta = rep(0, n.age-1)
  )
  
  Multinomial.obj <- MakeADFun(tmb.data,parameters,random=c("theta"), random.start=expression(last.par[random]),
                               DLL="Multinomial",inner.control=list(maxit=500,trace=F),silent = TRUE)  
  
  Multinomial.opt<-nlminb(Multinomial.obj$par, Multinomial.obj$fn, Multinomial.obj$gr,
                          control = list(trace=0,iter.max=10000,
                                         eval.max=10000,
                                         sing.tol=1e-20))
  
  Multinomial.rep = Multinomial.obj$report()
  
  sim.multinomial[iter,] = Multinomial.rep$p  
  
  ########################################################################
  parameters <- list( 
    log_phi = log(10),
    logit_ar = qlogis(0.5),
    log_sigmaP=log(0.1),
    theta = rep(0, n.age-1)
  )
  
  Diri_multinomial.obj <- MakeADFun(tmb.data,parameters,random=c("theta"), random.start=expression(last.par[random]),
                                    DLL="Multi_Diri",inner.control=list(maxit=500,trace=F),silent = TRUE)  
  
  Diri_multinomial.opt<-nlminb(Diri_multinomial.obj$par,Diri_multinomial.obj$fn, Diri_multinomial.obj$gr,
                               control = list(trace=0,iter.max=10000,
                                              eval.max=10000,
                                              sing.tol=1e-20))
  
  Diri_multinomial.rep = Diri_multinomial.obj$report()
  
  sim.Diri_multinomial[iter,] = Diri_multinomial.rep$p
  
  #################################################################################
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    last_not0 = last_not0-1,
    Xop=Xo.guess
  )
  
  parameters <- list( 
    #log_sigmaP=log(0.2),
    log_lambda0 = 2,
    logit_ar = rep(0,1), 
    log_sd = log(0.2),
    logit_ar_data = 0, 
    log_sd_data = log(0.2),
    theta = rep(0,n.age-1)
  )
  
  error_r=try( ( delta_clr_fix_datacorr.obj <- MakeADFun(tmb.data,parameters,random=c("theta"),
                                                         random.start=expression(last.par[random]), 
                                                         DLL="delta_clr_fix_datacorr", 
                                                         inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_clr_fix_datacorr.obj$env$tracemgc <- FALSE
  
  error_r=try( ( delta_clr_fix_datacorr.opt<-nlminb(delta_clr_fix_datacorr.obj$par,delta_clr_fix_datacorr.obj$fn,
                                                    delta_clr_fix_datacorr.obj$gr, 
                                                    control = list(trace=0,iter.max=10000,
                                                                   eval.max=10000,
                                                                   sing.tol=1e-20)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_clr_fix_datacorr.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_clr_fix_datacorr.obj$gr(delta_clr_fix_datacorr.opt$par)
  
  delta_clr_fix_datacorr.rep = delta_clr_fix_datacorr.obj$report() 
  # delta_clr_fix_datacorr.sdrep = sdreport(delta_clr_fix_datacorr.obj)
  # summary(delta_clr_fix_datacorr.sdrep)
  
  p_true_save_clr[iter,] = delta_clr_fix_datacorr.rep$p_true
  
  ###################################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    Po_ratio = Po_ratio,
    last_not0 = last_not0-1,
    Xop=Xo.guess
  )
  
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 2,
    logit_ar = rep(0,1), 
    log_sd = log(0.2),
    logit_ar_data = 0, 
    log_sd_data = log(0.2),
    theta = rep(0,n.age-1)
  )
  
  error_r=try( ( delta_mL_datacorr.obj <- MakeADFun(tmb.data,parameters,random=c("theta"),
                                                    random.start=expression(last.par[random]), 
                                                    DLL="delta_mL_datacorr", 
                                                    inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_mL_datacorr.obj$env$tracemgc <- FALSE
  
  error_r=try( ( delta_mL_datacorr.opt<-nlminb(delta_mL_datacorr.obj$par,delta_mL_datacorr.obj$fn,
                                               delta_mL_datacorr.obj$gr, 
                                               control = list(trace=0,iter.max=10000,
                                                              eval.max=10000,
                                                              sing.tol=1e-20)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_mL_datacorr.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_mL_datacorr.obj$gr(delta_mL_datacorr.opt$par)
  
  delta_mL_datacorr.rep = delta_mL_datacorr.obj$report() 
  # delta_mL_datacorr.sdrep = sdreport(delta_mL_datacorr.obj)
  # summary(delta_mL_datacorr.sdrep)
  
  p_true_save_ml[iter,] = delta_mL_datacorr.rep$p_true
  ##########################################################################
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 2,
    logit_ar0 = rep(0,1), 
    log_sd0 = log(0.2),
    logit_ar = rep(0,1), 
    log_sd = log(0.2),
    logit_ar_data = 0, 
    log_sd_data = log(0.2),
    p_mean = 0,
    theta = rep(0,n.age-1),
    dev_p0 = rep(0,n.age)
  )
  
  error_r=try( ( delta_mL_datacorr_logitp0.obj <- MakeADFun(tmb.data,parameters,random=c("theta","dev_p0"),
                                                            random.start=expression(last.par[random]), 
                                                            DLL="delta_mL_datacorr_logitp0", 
                                                            inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_mL_datacorr_logitp0.obj$env$tracemgc <- FALSE
  
  error_r=try( ( delta_mL_datacorr_logitp0.opt<-nlminb(delta_mL_datacorr_logitp0.obj$par,delta_mL_datacorr_logitp0.obj$fn,
                                                       delta_mL_datacorr_logitp0.obj$gr, 
                                                       control = list(trace=0,iter.max=10000,
                                                                      eval.max=10000,
                                                                      sing.tol=1e-20)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_mL_datacorr_logitp0.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_mL_datacorr_logitp0.obj$gr(delta_mL_datacorr_logitp0.opt$par)
  
  delta_mL_datacorr_logitp0.rep = delta_mL_datacorr_logitp0.obj$report() 
  # delta_mL_datacorr_logitp0.sdrep = sdreport(delta_mL_datacorr_logitp0.obj)
  # summary(delta_mL_datacorr_logitp0.sdrep)
  
  p_true_save_ml_logitp0[iter,] = delta_mL_datacorr_logitp0.rep$p_true
  
  #######################################################################################
  X = Po
  Y = sqrt(X)
  theta = matrix(NA,nrow=n.year,ncol=(n.age-1))
  
  #### Calculating theta by transforming X -> Y -> theta
  theta[,1] = acos(Y[,1])   ### i=1
  for (j in 1:n.year){
    for (i in 2:(n.age-1)){
      temp = Y[j,i]/prod(sin(theta[j,1:(i-1)]))
      if(temp>=1){temp = 0.999999}
      theta[j,i] = acos(temp)
      # print(i)
    }
  }
  
  
  tmb.data = list(
    theta=theta
  )
  
  
  parameters <- list( 
    log_sigmaP1=log(0.2),                  
    logit_ar1 = log(0.9/0.1),
    log_sigmaP2=log(0.2),                  
    logit_ar2 = log(0.9/0.1), 
    mu = rep(0,n.age-1)
  )
  
  DRHT.obj <- MakeADFun(tmb.data,parameters,random=c("mu"), 
                        DLL="DRHT",inner.control=list(maxit=500,trace=F),silent = TRUE)  
  
  DRHT.opt<-nlminb(DRHT.obj $par, DRHT.obj$fn, , DRHT.obj $gr,
                   control = list(trace=0,eval.max=2000,iter.max=1000))
  
  DRHT.rep = DRHT.obj$report() 
  
  theta_hat = exp(DRHT.rep$mu)
  Y_hat =rep(0, n.age)
  Y_hat[1]= cos(theta_hat[1])
  
  for (j in 2:(n.age-1)){
    Y_hat[j] = prod(sin(theta_hat[1:(j-1)]))*cos(theta_hat[j])
  }
  Y_hat[n.age] = prod(prod(sin(theta_hat[1:(n.age-1)])))
  X_hat = Y_hat^2
  
  p_true_save_DRHT[iter,] = X_hat
  
  ###################################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    Po_ratio = Po_ratio,
    last_not0 = last_not0-1,
    Xop=Xo.guess
  )
  
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 2,
    logit_ar = 0, 
    log_sd = log(0.2),
    log_phi = 0,
    theta = rep(0,n.age-1)
  )
  
  error_r=try( ( delta_Dirichlet.obj <- MakeADFun(tmb.data,parameters,random=c("theta"),
                                                  random.start=expression(last.par[random]), 
                                                  DLL="delta_Dirichlet", 
                                                  inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_Dirichlet.obj$env$tracemgc <- FALSE
  
  error_r=try( ( delta_Dirichlet.opt<-nlminb(delta_Dirichlet.obj$par,delta_Dirichlet.obj$fn,
                                             delta_Dirichlet.obj$gr, 
                                             control = list(trace=0,iter.max=10000,
                                                            eval.max=10000,
                                                            sing.tol=1e-20)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_Dirichlet.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_Dirichlet.obj$gr(delta_Dirichlet.opt$par)
  
  delta_Dirichlet.rep = delta_Dirichlet.obj$report() 
  # delta_Dirichlet.sdrep = sdreport(delta_Dirichlet.obj)
  # summary(delta_Dirichlet.sdrep)
  
  p_true_save_Dirichlet[iter,] = delta_Dirichlet.rep$p_true
  
  ###################################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    Po_ratio = Po_ratio,
    last_not0 = last_not0-1,
    Xop=Xo.guess
  )
  
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 2,
    logit_ar0 = rep(0,1), 
    log_sd0 = log(0.2),
    logit_ar = 0, 
    log_sd = log(0.2),
    log_phi = 0,
    theta = rep(0,n.age-1),
    dev_p0 = rep(0,n.age)
  )
  
  error_r=try( ( delta_Dirichlet_AR1_0.obj <- MakeADFun(tmb.data,parameters,random=c("theta","dev_p0"),
                                                        random.start=expression(last.par[random]), 
                                                        DLL="delta_Dirichlet_AR1_0", 
                                                        inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_Dirichlet_AR1_0.obj$env$tracemgc <- FALSE
  
  error_r=try( ( delta_Dirichlet_AR1_0.opt<-nlminb(delta_Dirichlet_AR1_0.obj$par,delta_Dirichlet_AR1_0.obj$fn,
                                                   delta_Dirichlet_AR1_0.obj$gr, 
                                                   control = list(trace=0,iter.max=10000,
                                                                  eval.max=10000,
                                                                  sing.tol=1e-20)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_Dirichlet_AR1_0.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_Dirichlet_AR1_0.obj$gr(delta_Dirichlet_AR1_0.opt$par)
  
  delta_Dirichlet_AR1_0.rep = delta_Dirichlet_AR1_0.obj$report() 
  # delta_Dirichlet_AR1_0.sdrep = sdreport(delta_Dirichlet_AR1_0.obj)
  # summary(delta_Dirichlet_AR1_0.sdrep)
  
  p_true_save_Dirichlet_0AR1[iter,] = delta_Dirichlet_AR1_0.rep$p_true
  
  ####################################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  Zo = matrix(NA,nrow=dim(Po)[1],ncol=dim(Po)[2]-1)
  for (i in 1:dim(Zo)[1]) {
    s_now = 0
    for (j in 1:dim(Zo)[2]) {
      if(s_now<1){
        Zo[i,j] = Po[i,j]/(1-s_now)
        s_now = s_now + Po[i,j]
      }
    }
  }
  
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    Po_ratio = Po_ratio,
    last_not0 = last_not0-1,
    Xop=Xo.guess,
    Zo = Zo
  )
  
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 0,
    logit_ar0 = rep(0,1), 
    log_sd0 = log(0.2),
    logit_ar = 0, 
    log_sd = log(0.2),
    # log_phia = log(3),
    logit_ar_b = 0, 
    log_sd_b = log(0.2),
    # log_phib = log(3),
    lai = rep(0,n.age-1),
    lbi = rep(0,n.age-1),
    dev_p0 = rep(0,n.age)
  )
  
  map = list(
    log_phia=as.factor(c(NA)),
    log_phib=as.factor(c(NA))
    # log_phia=as.factor(c("a")),
    # log_phib=as.factor(c("a"))
    # lbi = factor(rep(NA,n.age-2))
  )
  
  error_r=try( ( delta_generalized_Dirichlet_AR1_0_zi.obj <- MakeADFun(tmb.data,parameters,
                                                                       random=c("lai","lbi","dev_p0"),
                                                                       #map=map,
                                                                       random.start=expression(last.par[random]), 
                                                                       DLL="delta_generalized_Dirichlet_AR1_0_zi", 
                                                                       inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_generalized_Dirichlet_AR1_0_zi.obj$env$tracemgc <- FALSE
  
  # delta_generalized_Dirichlet_AR1_0_zi.obj$fn()
  # delta_generalized_Dirichlet_AR1_0_zi.obj$gr()
  
  error_r=try( ( delta_generalized_Dirichlet_AR1_0_zi.opt<-nlminb(delta_generalized_Dirichlet_AR1_0_zi.obj$par,
                                                                  delta_generalized_Dirichlet_AR1_0_zi.obj$fn,
                                                                  delta_generalized_Dirichlet_AR1_0_zi.obj$gr, 
                                                                  control = list(trace=0,iter.max=10000,
                                                                                 eval.max=10000,
                                                                                 sing.tol=1e-20)) ), 
               silent=TRUE )
  
  # delta_generalized_Dirichlet_AR1_0_zi.obj$gr(delta_generalized_Dirichlet_AR1_0_zi.opt$par)
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_generalized_Dirichlet_AR1_0_zi.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_generalized_Dirichlet_AR1_0_zi.obj$gr(delta_generalized_Dirichlet_AR1_0_zi.opt$par)
  
  delta_generalized_Dirichlet_AR1_0_zi.rep = delta_generalized_Dirichlet_AR1_0_zi.obj$report() 
  # delta_Dirichlet.sdrep = sdreport(delta_generalized_Dirichlet_AR1_0_zi.obj)
  # summary(delta_Dirichlet.sdrep)
  
  p_true_save_generalized_Dirichlet_AR1_0[iter,] = delta_generalized_Dirichlet_AR1_0_zi.rep$p_true
  
  ####################################################################################
  Xo.guess = Xop 
  #Xo.guess = rep(200,n.year)
  
  tmb.data = list(
    Po=Po,
    count_not0 = count_not0,
    Po_ratio = Po_ratio,
    last_not0 = last_not0-1,
    Xop=Xo.guess,
    Zo = Zo
  )
  
  parameters <- list( 
    # log_sigmaP=log(0.2),
    log_lambda0 = 0,
    logit_ar = 0, 
    log_sd = log(0.2),
    # log_phia = log(3),
    logit_ar_b = 0, 
    log_sd_b = log(0.2),
    # log_phib = log(3),
    lai = rep(0,n.age-1),
    lbi = rep(0,n.age-1)
  )
  
  map = list(
    log_phia=as.factor(c(NA)),
    log_phib=as.factor(c(NA))
    # log_phia=as.factor(c("a")),
    # log_phib=as.factor(c("a"))
    # lbi = factor(rep(NA,n.age-2))
  )
  
  error_r=try( ( delta_generalized_Dirichlet.obj <- MakeADFun(tmb.data,parameters,
                                                              random=c("lai","lbi"),
                                                              #map=map,
                                                              random.start=expression(last.par[random]), 
                                                              DLL="delta_generalized_Dirichlet", 
                                                              inner.control=list(maxit=10000,trace=F)) ), 
               silent=TRUE )
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  delta_generalized_Dirichlet.obj$env$tracemgc <- FALSE
  
  # delta_generalized_Dirichlet.obj$fn()
  # delta_generalized_Dirichlet.obj$gr()
  
  error_r=try( ( delta_generalized_Dirichlet.opt<-nlminb(delta_generalized_Dirichlet.obj$par,
                                                         delta_generalized_Dirichlet.obj$fn,
                                                         delta_generalized_Dirichlet.obj$gr, 
                                                         control = list(trace=0,iter.max=10000,
                                                                        eval.max=10000,
                                                                        sing.tol=1e-20)) ), 
               silent=TRUE )
  
  # delta_generalized_Dirichlet.obj$gr(delta_generalized_Dirichlet.opt$par)
  
  if ('try-error' %in% class(error_r)) {
    error_count = error_count + 1
    next
  }
  
  if (max(abs(delta_generalized_Dirichlet.obj$gr()))>0.01) {
    error_count = error_count + 1
    next
  }
  
  # delta_generalized_Dirichlet.obj$gr(delta_generalized_Dirichlet.opt$par)
  
  delta_generalized_Dirichlet.rep = delta_generalized_Dirichlet.obj$report() 
  # delta_Dirichlet.sdrep = sdreport(delta_generalized_Dirichlet.obj)
  # summary(delta_Dirichlet.sdrep)
  
  p_true_save_generalized_Dirichlet[iter,] = delta_generalized_Dirichlet.rep$p_true
  
  ####################################################################################
  
  print(iter)
  iter = iter+1
  
}










