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

parameterSummary = summary(T_modelTR, c("R0","r_lock_1","m_lock_1","shift_lock_1","gamma_s","gamma_H","gamma_ICU","eps_H_ICU",
                                        "eps_H_x","eps_ICU_x","r_d_s","r_d_a"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975))
write.table(parameterSummary$summary, file = paste0(path2save,"/CSVS/parameters_fit_mrelax_",m_relax_in,".csv"))

pp3 = c("predicted_daily_cases","predicted_daily_deaths","predicted_current_icu")
model_output_all = summary(T_modelTR,pp3)[[1]] %>%
  tbl_df() %>%
  mutate(t=rep(1:data_list$S,3),
         date=date_data+t-1,
         eta="-100%",
         type=rep(pp3,each=data_list$S)) %>%
  mutate(populasyon=factor(type,levels=pp3,
                           labels=c("Gunluk Vaka Sayisi","Gunluk Vefat Sayisi","Yogun Bakim Hasta Sayisi")))

write.table(model_output_all, file = paste0(path2save,"/CSVS/output_all_mrelax_",m_relax_in,".csv"))

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
ggsave(paste0(path2save,"/FIGS/figure_all_mrelax_",m_relax_in,".png"),width = 10, height = 5)

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

