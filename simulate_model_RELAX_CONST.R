# Setup
rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# path2save = paste0("OUT_",format(Sys.time(), "%d_%b_%Y"));
# dir.create(path2save)
# dir.create(paste0(path2save,"/CSVS/"))
# dir.create(paste0(path2save,"/FIGS/"))
# dir.create(paste0(path2save,"/RDATA/"))
library(matlab)
library(deSolve)

load("/Users/burcu/Dropbox/corona-tr-modeling/OUT_10_Jun_2020/RDATA/T_modelTR_mrelax_4.RData")
parameterSummary = summary(T_modelTR, c("R0","tau","gamma_s","gamma_H","gamma_ICU","eps_H_ICU",
                                        "eps_H_x","eps_ICU_x","r_d_s","r_d_a","shift_lock_1","m_lock_1",
                                        "r_lock_1","pi"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))




source("prepare_data.R")
K   = test_fit_vec_1[1];
mu  = test_fit_vec_1[2];
sig = test_fit_vec_1[3];

############

params_a = as.numeric(unlist(parameterSummary$summary[1:13,6]))
pi_med   = as.numeric(unlist(parameterSummary$summary[14,6]))
popsize  = 84180493;

covid.model <- function (t, y, params) {
  
  S   = y[1];
  E   = y[2];
  I   = y[3];
  H   = y[4];
  ICU = y[5];
  R   = y[6];
  X   = y[7];
  C   = y[8];
  
  R0           = params[1]
  tau          = params[2]
  gamma_s      = params[3]
  gamma_H      = params[4]
  gamma_ICU    = params[5]
  eps_H_ICU    = params[6]
  eps_H_x      = params[7]
  eps_ICU_x    = params[8]
  r_d_s        = params[9]
  r_d_a        = params[10]
  
  shift_lock   = params[11]
  m_lock       = params[12]
  r_lock       = params[13]
  t_lock       = params[14]
  
  shift_relax  = params[15]
  m_relax      = params[16]
  t_normal     = params[17]
  r_end        = params[18]
  
  K   = params[19]
  mu  = params[20]
  sig = params[21]

  piVal    = 3.14;
  
  p_lock   = r_lock+(1-r_lock)/(1+exp(m_lock*(t-t_lock-shift_lock)));
  p_relax  = r_lock+1./(1/(r_end-r_lock)+exp(-m_relax*(t-t_normal-shift_relax)));
  p_asymp  = K*(1/(t*sig*sqrt(2*piVal)))*exp(-((log(t)-mu)/(sqrt(2)*sig))^2);
  
  if(is.nan(p_asymp)){
    p_asymp = 0;
  }
  if(p_asymp<0){
    p_asymp = 0;
  }
  
  if(p_relax>p_lock){
    coeffR = p_relax;
  }
  else{
    coeffR = p_lock;
  }
  
  print(c(t,coeffR))
  
  dSdt = - coeffR*S*(R0*gamma_s*I); 
  dEdt = + coeffR*S*(R0*gamma_s*I) - tau*E; 
  dIdt = + tau*E - gamma_s*I;
  dHdt = + r_d_s*gamma_s*I - gamma_H*H; 
  dICUdt = + gamma_H*eps_H_ICU*H - gamma_ICU*ICU;
  dRdt = + gamma_H*(1-eps_H_ICU-eps_H_x)*H + gamma_ICU*(1-eps_ICU_x)*ICU + (1-r_d_s)*gamma_s*I; 
  dXdt = + gamma_H*eps_H_x*H + gamma_ICU*eps_ICU_x*ICU; 
  dCdt = + (r_d_s + p_asymp*r_d_a)*gamma_s*I;
  
  dydt = c(dSdt,dEdt,dIdt,dHdt,dICUdt,dRdt,dXdt,dCdt)
  list(dydt)
}

ystart <- c(S=1E-9+(1-pi_med),
            E=1E-9+(pi_med),
            I=1E-9,
            H=1E-9,
            ICU=1E-9,
            R=1E-9,
            X=1E-9,
            C=1E-9) # initial conditions

times       =  seq(from=0,to=180,by=1) # returns a sequence
t_lock      = 9;
t_normal   = 360;

# optimistic
shift_relax= -200;
m_relax    = 0.05;
r_end      = 0.25;

params_b = c(t_lock,shift_relax,m_relax,t_normal,r_end)
params   = c(params_a,params_b,K,mu,sig)

out_o <- as.data.frame(
  ode(
    func=covid.model,
    y=ystart,
    times=times,
    parms=params
  )
)

out_o[,2:9] = popsize*out_o[,2:9];
dC    = c(0,diff(out_o$C));
dX    = c(0,diff(out_o$X));
out_o = tibble(out_o,dC,dX)

# pessimistic
shift_relax= -240;
m_relax    = 0.1;
r_end      = 0.25;

params_b = c(t_lock,shift_relax,m_relax,t_normal,r_end)
params   = c(params_a,params_b,K,mu,sig)

out_p <- as.data.frame(
  ode(
    func=covid.model,
    y=ystart,
    times=times,
    parms=params
  )
)


out_p[,2:9] = popsize*out_p[,2:9];
dC    = c(0,diff(out_p$C));
dX    = c(0,diff(out_p$X));
out_p = tibble(out_p,dC,dX)

par(mfrow=c(2,2))
plot(dC~time,data=out_o,type='l',main="C",col="blue")
lines(dC~time,data=out_p,type='l',main="C",col="blue")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
points(daily_cases_data)
plot(dX~time,data=out_o,type='l',main="X",col="red")
lines(dX~time,data=out_p,type='l',main="X",col="red")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
points(daily_deaths_data)
plot(ICU~time,data=out_o,type='l',main="ICU",col="pink")
lines(ICU~time,data=out_p,type='l',main="ICU",col="pink")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
points(icu_data)



