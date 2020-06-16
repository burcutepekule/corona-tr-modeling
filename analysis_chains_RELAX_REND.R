library("bayesplot")
library("ggplot2")
posterior  = as.array(T_modelTR)

color_scheme_set("blue")
mcmc_res=mcmc_dens(posterior, pars = c("R0","r_lock_1","r_end","m_lock_1","m_relax_1",
                                       "gamma_s","gamma_H","gamma_ICU","eps_H_ICU",
                              "eps_H_x","eps_ICU_x","r_d_s","r_d_a"),
          facet_args = list(ncol = 4, strip.position = "left"))
ggsave(paste0(path2save,"/FIGS/posteriors.png"),width = 10, height = 6)

color_scheme_set("viridis")
mcmc_vals=mcmc_trace(posterior,pars = c("R0","r_lock_1","r_end","m_lock_1","m_relax_1",
                                        "gamma_s","gamma_H","gamma_ICU","eps_H_ICU",
                              "eps_H_x","eps_ICU_x","r_d_s","r_d_a"),
           facet_args = list(ncol = 3, strip.position = "left"))
ggsave(paste0(path2save,"/FIGS/chains.png"),width = 10, height = 6)