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

m_relax_in = 0 # slope of relaxation -> if 0, no relaxation

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
  K   = test_fit_vec_1[1],
  mu  = test_fit_vec_1[2],
  sig = test_fit_vec_1[3],
  D=as.numeric(date_end-date_data+1),
  k_daily_cases  = daily_cases_data,
  k_icu          = icu_data,
  k_daily_deaths = daily_deaths_data,
  
  p_pi        = c(1,999),
  p_R0        = c(3,1),
  p_tau       = c(2,3.2),
  p_gamma_s   = c(2,3.2),
  p_gamma_H   = c(2,3),
  p_gamma_ICU = c(2,3),
  p_eps_H_ICU = c(2,3),
  p_eps_H_x   = c(2,10),
  p_eps_ICU_x = c(2,10),
  p_r_d_s     = c(2,10),
  p_r_d_a     = c(2,10),
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

load("/Users/burcu/Dropbox/corona-tr-modeling/OUT_31_May_2020/RDATA/T_modelTR_mrelax_0.RData")
T_modelTR_old = T_modelTR
model_end = ymd("2020-05-31")
D_old=as.numeric(model_end-date_data+1)+days2add

pp = c("predicted_daily_cases","predicted_daily_deaths","predicted_current_icu")
model_output_fitted_old = summary(T_modelTR_old,pp)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:D_old,3),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp,each=D_old)) %>%
  mutate(populasyon=factor(type,levels=pp,
                           labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi")))

D_data = as.numeric(date_end-date_data+1)
populasyon_data=factor(rep(pp,each=D_data),levels=pp,
                  labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi"))

agg_data         = tibble(allDates_agg,daily_cases_data,daily_deaths_data,icu_data);
observed_data    = c(daily_cases_data,daily_deaths_data,icu_data);
observed_dates   = c(allDates_agg,allDates_agg,allDates_agg);
agg_data_all     = tibble(all_dates=observed_dates,all_data=observed_data,populasyon=populasyon_data);


old_data_length     = as.numeric(model_end-date_data+1);
populasyon_data_old =factor(rep(pp,each=old_data_length),levels=pp,
                       labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi"))
observed_data_old   = c(daily_cases_data[1:old_data_length],daily_deaths_data[1:old_data_length],icu_data[1:old_data_length]);
observed_dates_old  = c(allDates_agg[1:old_data_length],allDates_agg[1:old_data_length],allDates_agg[1:old_data_length]);
agg_data_all_old    = tibble(all_dates=observed_dates_old,all_data=observed_data_old,populasyon=populasyon_data_old);

# # 
ggplot() +
  geom_ribbon(data=model_output_fitted_old,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populasyon),alpha=.5) +
  geom_line(data=model_output_fitted_old,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populasyon ,scales="free",nrow=2) +
  geom_vline(xintercept=date_control_3,linetype=2) +
  geom_vline(xintercept=model_end,linetype=1) +
  geom_vline(xintercept=date_end,linetype=1) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_point(data=agg_data_all,aes(x=all_dates,y=all_data,fill=populasyon),shape=21,size=1,colour = "black", fill = "black") +
  geom_point(data=agg_data_all_old,aes(x=all_dates,y=all_data,fill=populasyon),shape=21,size=1,colour = "black", fill = "white") +
  labs(x="Gun",y="Birey Sayisi")
# scale_y_continuous(trans = 'log10')
ggsave(paste0(path2save,"/FIGS/figure_fit_mrelax_",m_relax_in,"_RETRO.png"),width = 10, height = 5)
