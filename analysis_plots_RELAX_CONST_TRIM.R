# Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path2save = paste0("OUT_",format(Sys.time(), "%d_%b_%Y"));
dir.create(path2save)
dir.create(paste0(path2save,"/CSVS/"))
dir.create(paste0(path2save,"/FIGS/"))
dir.create(paste0(path2save,"/RDATA/"))


library(rstan)
library(zoo)
library(Rcpp)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("prepare_data.R")
source("setup.R")

days2add     = 240; #ADDITIONAL DAYS FOR SIMULATION
date_simul   = date_end + days2add;

# SCENARIO 1
m_relax_in   = 5; #indicator for using relaxation, change this index given scenario
date_normal  = ymd("2020-08-31");
t_normal     = as.numeric(date_normal-date_data+1);
t_relax      = as.numeric(date_relax-date_data+1);
t_mid        = round(0.5*(t_normal-t_relax));
srelax       = (t_normal-t_relax);
mrelax       = 0.1;

data_list = list(
  pop_t    = pop_size,
  tlock_1  = as.numeric(date_control_1-date_data+1),
  trelax_1 = t_mid,
  mrelax_1 = mrelax,
  srelax_1 = srelax, 
  r_end    = 1, #back to complete normal
  K   = test_fit_vec_1[1],
  mu  = test_fit_vec_1[2],
  sig = test_fit_vec_1[3],
  mult= 0, #dummy control
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
  ts_pred=1:as.numeric(date_simul-date_data+1)
  
)

########################
load(paste0("/Users/burcu/Dropbox/corona-tr-modeling/OUT_10_Jun_2020/RDATA/T_modelTR_mrelax_",m_relax_in,".RData"))

pp = c("fitted_k_daily_cases","fitted_k_daily_deaths","fitted_k_icu")
model_output_fitted = summary(T_modelTR,pp)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$D,3),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp,each=data_list$D)) %>%
  mutate(populasyon=factor(type,levels=pp,
                           labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi")))

agg_data         = tibble(allDates_agg,daily_cases_data,daily_deaths_data,icu_data);
observed_data    = c(daily_cases_data,daily_deaths_data,icu_data);
observed_dates   = c(allDates_agg,allDates_agg,allDates_agg);
agg_data_all     = tibble(all_dates=observed_dates,all_data=observed_data,populasyon=model_output_fitted$populasyon);
# # 
ggplot() +
  geom_ribbon(data=model_output_fitted,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populasyon),alpha=.5) +
  geom_line(data=model_output_fitted,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populasyon ,scales="free",nrow=2) +
  geom_vline(xintercept=date_control_3,linetype=2) +
  geom_vline(xintercept=date_relax,linetype=2) +
  geom_vline(xintercept=date_end,linetype=3) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_point(data=agg_data_all,aes(x=all_dates,y=all_data,fill=populasyon),shape=21,size=2,colour = "black", fill = "white") +
  labs(x="Gun",y="Birey Sayisi")
# scale_y_continuous(trans = 'log10')
ggsave(paste0(path2save,"/FIGS/figure_fit_mrelax_",m_relax_in,".png"),width = 10, height = 5)

# parameterSummary = summary(T_modelTR, c("R0","r_lock_1","m_lock_1","shift_lock_1","gamma_s","gamma_H","gamma_ICU","eps_H_ICU",
#                                         "eps_H_x","eps_ICU_x","r_d_s","r_d_a"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
# write.table(parameterSummary$summary, file = paste0(path2save,"/CSVS/parameters_fit_mrelax_",m_relax_in,".csv"))

pp3 = c("predicted_daily_cases","predicted_daily_deaths","predicted_current_icu")
model_output_all = summary(T_modelTR,pp3)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,3),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp3,each=data_list$S)) %>%
  mutate(populasyon=factor(type,levels=pp3,
                           labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi")))

write.table(model_output_all, file = paste0(path2save,"/CSVS/output_all_short_mrelax_",m_relax_in,".csv"))

date_plot_until  = ymd("2020-08-01");
plot_length      = as.numeric(date_plot_until-date_data+1);

agg_data         = tibble(allDates_agg,daily_cases_data,daily_deaths_data,icu_data);
p_cases          = model_output_all %>% filter(type=="predicted_daily_cases")
p_deaths         = model_output_all %>% filter(type=="predicted_daily_deaths")
p_icu            = model_output_all %>% filter(type=="predicted_current_icu")

p_cases = p_cases[1:plot_length,]
p_deaths= p_deaths[1:plot_length,]
p_icu   = p_icu[1:plot_length,]

model_output_all = bind_rows(p_cases,p_deaths,p_icu)

ggplot() +
  geom_ribbon(data=model_output_all,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=populasyon),alpha=.5) +
  # geom_ribbon(data=model_output_all,aes(x=date,ymin=`25%`,ymax=`75%`,fill=populasyon),alpha=.5) +
  geom_line(data=model_output_all,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ populasyon ,scales="free",nrow=2) +
  geom_vline(xintercept=date_control_3,linetype=2) +
  geom_vline(xintercept=date_end,linetype=2) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_point(data=agg_data_all,aes(x=all_dates,y=all_data,fill=populasyon),shape=21,size=2,colour = "black", fill = "white") +
  labs(x="Gun",y="Birey Sayisi")
# scale_y_continuous(trans = 'log10')
ggsave(paste0(path2save,"/FIGS/figure_all_short_mrelax_",m_relax_in,".png"),width = 10, height = 5)


##
pp5 = c("RE_raw")
RE_output_all = summary(T_modelTR,pp5)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,1),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp5,each=data_list$S)) %>%
  mutate(Re=factor(type,levels=pp5,
                   labels=c("RE")))

selectedRows_RE <- RE_output_all[grep("RE", RE_output_all$Re), ]
selectedRows_RE = 

ggplot() +
  geom_ribbon(data=selectedRows_RE,aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=Re),alpha=.5) +
  geom_line(data=selectedRows_RE,aes(x=date,y=`50%`),colour="black") +
  facet_wrap(~ Re ,scales="free",nrow=1) +
  geom_vline(xintercept=date_control_3,linetype=2) +
  geom_vline(xintercept=date_end,linetype=2) +
  scale_colour_manual(values=c("grey20","black"),guide=FALSE) +
  scale_alpha_manual(values=c(1,0),guide=FALSE) +
  scale_x_date(date_breaks="2 weeks",date_labels="%b %d") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x="Gun",y="Re")
# scale_y_continuous(trans = 'log10')
ggsave(paste0(path2save,"/FIGS/figure_Re_mrelax_",m_relax_in,".png"),width = 8, height = 3)
