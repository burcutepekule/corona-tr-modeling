# Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path2save = paste0("OUT_",format(Sys.time(), "%d_%b_%Y"));
dir.create(path2save)
dir.create(paste0(path2save,"/CSVS/"))
dir.create(paste0(path2save,"/FIGS/"))
dir.create(paste0(path2save,"/RDATA/"))

args <- commandArgs(TRUE)
print(args)

m_relax_in = 0 # slope of relaxation

library(rstan)
library(zoo)
library(Rcpp)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("prepare_data.R")
source("setup.R")

days2add     = 45; #ADDITIONAL DAYS FOR SIMULATION
date_simul   = date_end + days2add;

data_list = list(
  pop_t=pop_size,
  tswitch_1=as.numeric(date_control_3-date_data+1),
  trelax=as.numeric(date_relax-date_data+1),
  r_end=1,
  m_relax=m_relax_in/100, 
  p1=test_fit_vec[1],
  p2=test_fit_vec[2],
  p3=test_fit_vec[3],
  p4=test_fit_vec[4],
  D=as.numeric(date_end-date_data+1),
  k_daily_cases  = daily_cases_data,
  k_icu          = icu_data,
  k_daily_deaths = daily_deaths_data,
  k_daily_recov  = recov_data,

  p_pi        = c(1,999),
  p_R0        = c(3,1),
  p_tau       = c(2,3.2),
  p_gamma_s   = c(2,3.2),
  p_gamma_H   = c(2,3),
  p_gamma_ICU = c(2,3),
  p_eps_H_ICU = c(2,3),
  p_eps_H_x   = c(2,10),
  p_eps_ICU_x = c(2,10),
  p_eps_H     = c(2,3),
  p_r_d_s     = c(2,10),
  p_r_d_a     = c(2,10),
  p_r_r_s     = c(1,1),
  p_r_r_a     = c(1,1),
  p_r_lock_1  = c(1,1),
  p_phi       = 1/100,

  t0=0,
  t_data=1,
  S=as.numeric(date_simul-date_data+1),
  E=days2add, 
  ts=1:as.numeric(date_end-date_data+1),
  ts_pred=1:as.numeric(date_simul-date_data+1),
  
  r_c=0 #dummy control
)
# # IF .rds not compiled (run in case of change in model)
M_model_TR     = stan_model("MODELS/model_TR.stan")
# M_model_TR = readRDS("MODELS/model_TR.rds")
####### DEBUG MODE
T_modelTR      = sampling(M_model_TR,data = data_list,iter=5,chains=1,init="random") 
####### FITTING
T_modelTR      = sampling(M_model_TR,data = data_list,warmup=150,iter=300,chains=8,init="random")
save(T_modelTR, file =paste0(path2save,"/RDATA/T_modelTR_mrelax_",m_relax_in,".RData"))

source("analysis_plots.R")
source("analysis_chains.R")


