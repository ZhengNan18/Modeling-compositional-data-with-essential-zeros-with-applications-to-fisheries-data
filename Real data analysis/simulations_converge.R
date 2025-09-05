library(TMB) 
library('ggplot2')  
library('gridExtra')
#library('ggcorrplot')

load("length_1996_2018_save")
head(length_1996_2018)
dim(length_1996_2018)
table(length_1996_2018[,144])
year_chosen = 2018
data_comp = length_1996_2018[length_1996_2018[,144]==year_chosen,]
data_comp = data_comp[,-(143:144)]
data_comp=data_comp[rowSums(data_comp!=0)>1,]
dim(data_comp)
data_proportion = colSums(length_1996_2018[length_1996_2018[,144]==year_chosen,-(143:144)])/
  sum(colSums(length_1996_2018[length_1996_2018[,144]==year_chosen,-(143:144)]))

n.age = dim(data_comp)[2] 
n.year=dim(data_comp)[1]  

# logitp0
compile("Multi_Diri.cpp")  

dyn.load(dynlib("Multi_Diri"))
#dyn.unload(dynlib("Multi_Diri"))

# logitp0
compile("delta_mL_datacorr_logitp0_jacob.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_mL_datacorr_logitp0_jacob"))
#dyn.unload(dynlib("delta_mL_datacorr_logitp0_jacob"))

####
compile("NBac_AR1.cpp")  

dyn.load(dynlib("NBac_AR1"))  
# dyn.load("NBac_AR1") 
#dyn.unload(dynlib("NBac_AR1")) 

#Dirichlet
compile("delta_Dirichlet_AR1_0.cpp",
        flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")

dyn.load(dynlib("delta_Dirichlet_AR1_0"))
#dyn.unload(dynlib("delta_Dirichlet_AR1_0"))

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
p_true_save_ml_logitp0_perturb = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_DRHT = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_Dirichlet = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_Dirichlet_0AR1 = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_AR1_0 = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_AR1_0_perturb = matrix(NA,nrow=nsim,ncol=n.age)

p_true_save_generalized_Dirichlet_perturb = matrix(NA,nrow=nsim,ncol=n.age)

set.seed(100)
################################################################################################
Xop = rowSums(data_comp) 
#Xop = rep(mean_tot_X,n.year)  

Xo = data_comp
Xr = Xo

Po = Xo/matrix(Xop,nrow=n.year,ncol=n.age,byrow=F)
rowSums(Po)

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
# Xo.guess = rep(100,n.year)

tmb.data = list(Po=Po,Xop=Xo.guess,upbound=0.9)

parameters_NBac_AR1 <- list( 
  log_sigmaP=log(0.1),                  
  logit_ar = 1, 
  log_sd_theta = log(0.2),
  logit_ar_theta = 0,
  theta = rep(0,n.age-1), 
  log_gam=matrix(0,nrow=n.year,ncol=n.age)
  # log_gamp=matrix(0,nrow=n.year,ncol=n.age)
  # log_gamt = rep(0, n.age)
)

NBac_AR1.obj <- MakeADFun(tmb.data,parameters_NBac_AR1,random=c("theta","log_gam"), 
                          DLL="NBac_AR1",inner.control=list(maxit=5000,trace=F),silent = TRUE)  

NBac_AR1.opt<-nlminb(NBac_AR1.obj$par,NBac_AR1.obj$fn,NBac_AR1.obj$gr,
                     control = list(trace=10,eval.max=2000,iter.max=1000))

NBac_AR1.rep = NBac_AR1.obj$report() 
# NBac_AR1.rep = NBac_AR1.obj$simulate()

# NBac_AR1.sdrep = sdreport(NBac_AR1.obj)

sim.nbac_AR1[iter,] = NBac_AR1.rep$p 
sim.nbac_AR1_true[iter,] = NBac_AR1.rep$p_true 

range(NBac_AR1.rep$gi_array)
NBac_AR1.rep$gi_array[,4]
# residual ##########################################################
Bin_qr = function(x,n,p) { 
  if(x==0){
    qr = runif(1, min = 0, max = pbinom(0, n, p))
  }else{
    qr = runif(1, min = pbinom(x-1, n, p), max = pbinom(x, n, p))
  }
}

res_vec=c()

set.seed(20)
for (i in 1:n.year) {
  data_now = round(Po[i,]*Xo.guess[i])
  total_now = sum(data_now)
  p_now = NBac_AR1.rep$rpm[i,]
  res_vec=c(res_vec,Bin_qr(data_now[1],total_now,p_now[1]))
  data_cum = data_now[1]
  prob_cum = p_now[1]
  for (j in 2:(n.age-1)) {
    res_vec=c(res_vec,Bin_qr(data_now[j],total_now-data_cum,p_now[j]/(1-prob_cum)))
    data_cum = data_cum + data_now[j]
    prob_cum = prob_cum + p_now[j]
  }
}

res_vec = qnorm(res_vec)
qqnorm(res_vec)
qqline(res_vec)

ks.test(res_vec,pnorm)

#################################################################################

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
  logit_ar = rep(0,1), 
  log_sd = log(0.2),
  logit_ar_data = 0, 
  log_sd_data = log(0.2),
  p_mean = 0,
  theta = rep(0,n.age-1),
  dev_p0 = rep(0,n.age)
)

map = list(
  p_mean = factor(NA)
)

error_r=try( ( delta_mL_datacorr_logitp0_jacob.obj <- MakeADFun(tmb.data,parameters,random=c("theta","dev_p0"),
                                                                map=map,
                                                                random.start=expression(last.par[random]), 
                                                                DLL="delta_mL_datacorr_logitp0_jacob", 
                                                                inner.control=list(maxit=10000,trace=F)) ), 
             silent=TRUE )

if ('try-error' %in% class(error_r)) {
  error_count = error_count + 1
  next
}

delta_mL_datacorr_logitp0_jacob.obj$env$tracemgc <- FALSE

error_r=try( ( delta_mL_datacorr_logitp0_jacob.opt<-nlminb(delta_mL_datacorr_logitp0_jacob.obj$par,
                                                           delta_mL_datacorr_logitp0_jacob.obj$fn,
                                                           delta_mL_datacorr_logitp0_jacob.obj$gr, 
                                                           control = list(trace=10,iter.max=10000,
                                                                          eval.max=10000,
                                                                          sing.tol=1e-20)) ), 
             silent=TRUE )

if ('try-error' %in% class(error_r)) {
  error_count = error_count + 1
  next
}

if (max(abs(delta_mL_datacorr_logitp0_jacob.obj$gr()))>0.01) {
  error_count = error_count + 1
  next
}

# delta_mL_datacorr_logitp0_jacob.obj$gr(delta_mL_datacorr_logitp0_jacob.opt$par)

delta_mL_datacorr_logitp0_jacob.rep = delta_mL_datacorr_logitp0_jacob.obj$report() 
# delta_mL_datacorr_logitp0_jacob.sdrep = sdreport(delta_mL_datacorr_logitp0_jacob.obj)
# summary(delta_mL_datacorr_logitp0_jacob.sdrep)

p_true_save_ml_logitp0[iter,] = delta_mL_datacorr_logitp0_jacob.rep$p_true

# residual with 0 ########################################
length(delta_mL_datacorr_logitp0_jacob.rep$p0)
p0_now = delta_mL_datacorr_logitp0_jacob.rep$p0
length(delta_mL_datacorr_logitp0_jacob.rep$theta)
theta_now = c(delta_mL_datacorr_logitp0_jacob.rep$theta,0)
rho_now = delta_mL_datacorr_logitp0_jacob.rep$ar_data
var_now = (1-rho_now^2)*delta_mL_datacorr_logitp0_jacob.rep$sd_data

# res_non0 = c()
set.seed(20)
res_vec = c()
for (i in 1:n.year) {
  non0_1 = TRUE
  for (j in 1:n.age) {
    if(Po_ratio[i,j]==0){
      res_vec = c(res_vec, runif(1, min = 0, max = p0_now[j]))
      if(!non0_1){
        dev_prev = rho_now*dev_prev
        var_current = var_current + var_now*(rho_now^(pow_now*2))
        pow_now = pow_now+1
      }
    }else{
      if(j==last_not0[i]){
        temp_now = runif(1, min = p0_now[j], max = 1)
        res_vec = c(res_vec, temp_now)
      }
      
      if(j<last_not0[i]){
        dev_current = log(Po_ratio[i,j])-theta_now[j]+theta_now[last_not0[i]]
      }
      if(non0_1){
        non0_1 = FALSE
        res_now = dev_current/sqrt(delta_mL_datacorr_logitp0_jacob.rep$sd_data)
      }else{
        if(j<last_not0[i]){
          res_now = (dev_current-rho_now*dev_prev)/sqrt(var_current)
        }
      }
      res_now = p0_now[j] + pnorm(res_now)*(1-p0_now[j])
      res_vec = c(res_vec,res_now)
      var_current = var_now
      pow_now = 1
      dev_prev = dev_current
    }
  }
}

res_vec = qnorm(res_vec)
qqnorm(res_vec)
qqline(res_vec)
ks.test(res_vec,pnorm)

################################################################################
Xo.guess = Xop 
#Xo.guess = rep(200,n.year)

tmb.data = list(Po=Po,Xop=Xo.guess)

parameters <- list( 
  log_phi = log(10),
  logit_ar = qlogis(0.5),
  log_sigmaP=log(0.1),
  theta = rep(0, n.age-1)
)

Diri_multinomial.obj <- MakeADFun(tmb.data,parameters,random=c("theta"), random.start=expression(last.par[random]),
                                  DLL="Multi_Diri",inner.control=list(maxit=500,trace=F),silent = TRUE)  

Diri_multinomial.opt<-nlminb(Diri_multinomial.obj$par,Diri_multinomial.obj$fn, Diri_multinomial.obj$gr,
                             control = list(trace=10,iter.max=10000,
                                            eval.max=10000,
                                            sing.tol=1e-20))

Diri_multinomial.rep = Diri_multinomial.obj$report()

# Diri_multinomial.sdrep = sdreport(Diri_multinomial.obj)
# summary(Diri_multinomial.sdrep)

sim.Diri_multinomial[iter,] = Diri_multinomial.rep$p
# residual ################################################
library(VGAM)

BB_qr = function(x,n,shape1,shape2) { 
  if(x==0){
    qr = runif(1, min = 0, max = pbetabinom.ab(0, n, shape1, shape2))
  }else{
    qr = runif(1, min = pbetabinom.ab(x-1, n, shape1, shape2), max = pbetabinom.ab(x, n, shape1, shape2))
  }
}

res_vec=c()
set.seed(20)
alpha_now = Diri_multinomial.rep$alpha
for (i in 1:n.year) {
  data_now = round(Po[i,]*Xo.guess[i])
  total_now = sum(data_now)
  alpha_cum = sum(alpha_now) - alpha_now[1]
  res_vec=c(res_vec,BB_qr(data_now[1],total_now,alpha_now[1],alpha_cum))
  data_cum = data_now[1]
  for (j in 2:(n.age-1)) {
    alpha_cum = alpha_cum - alpha_now[j]
    res_vec=c(res_vec,BB_qr(data_now[j],total_now-data_cum,alpha_now[j],alpha_cum))
    data_cum = data_cum + data_now[j]
  }
}

res_vec = qnorm(res_vec)
qqnorm(res_vec)
qqline(res_vec)

ks.test(res_vec,pnorm)

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
                                                 control = list(trace=10,iter.max=10000,
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

# residual with 0 ########################################
length(delta_Dirichlet_AR1_0.rep$p0)
p0_now = delta_Dirichlet_AR1_0.rep$p0
length(delta_Dirichlet_AR1_0.rep$alpha)
alpha_now = delta_Dirichlet_AR1_0.rep$alpha

# res_non0 = c()
set.seed(20)
res_vec = c()
for (i in 1:n.year) {
  non0_1 = TRUE
  cum_alpha = 0
  for (j in 1:n.age) {
    if(Po[i,j]!=0){
      cum_alpha = cum_alpha + alpha_now[j]
    }
  }
  for (j in 1:n.age) {
    if(Po[i,j]==0){
      res_vec = c(res_vec, runif(1, min = 0, max = p0_now[j]))
    }else{
      cum_alpha = cum_alpha - alpha_now[j]
      if(non0_1){
        non0_1 = FALSE
        cum_data = Po[i,j]
        res_vec = c( res_vec, p0_now[j] + pbeta(Po[i,j], alpha_now[j], cum_alpha)*(1-p0_now[j]) )
      }else{
        if(j==last_not0[i]){
          res_vec = c( res_vec, runif(1, min = p0_now[j], max = 1) )
        }else{
          res_vec = c( res_vec, p0_now[j] + pbeta(Po[i,j]/(1-cum_data), alpha_now[j], cum_alpha)*(1-p0_now[j]) )
          cum_data = cum_data + Po[i,j]
        }
      }
    }
  }
}
res_vec = qnorm(res_vec)
qqnorm(res_vec)
qqline(res_vec)
ks.test(res_vec,pnorm)

################################################################################################



