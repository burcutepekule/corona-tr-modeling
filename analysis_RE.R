# Setup
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path2save = paste0("OUT_",format(Sys.time(), "%d_%b_%Y"));


library(R0)
library(zoo)
library(Rcpp)
library(EpiEstim)
library(ggplot2)
library(incidence)
source("prepare_data.R")
source("setup.R")


plot(as.incidence(daily_cases_data[-1], dates = allDates_agg[-1]))

m_relax_in = 1
# load(paste0(path2save,"/RDATA/T_modelTR_mrelax_",m_relax_in,".RData"))
load("/Users/burcu/Dropbox/corona-tr-modeling/OUT_24_Jun_2020/RDATA/T_modelTR_mrelax_2.RData")

     
parameterSummary = summary(T_modelTR, c("gamma_s","tau"), probs = c(0.025, 0.25, 0.50, 0.75, 0.975));

list_of_draws <- rstan::extract(T_modelTR)
si=(1/(list_of_draws$gamma_s)+1/(list_of_draws$tau))

res_parametric_si <- estimate_R(daily_cases_data[-1], 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = mean(si), 
                                  std_si = 10*sd(si)))
)

head(res_parametric_si$R)
res_parametric_si$dates=allDates_agg[-1];
png(filename=paste0(path2save,"/FIGS/RE_estimate.png"),
    width     = 6.25,
    height    = 6.25,
    units     = "in",
    res       = 1200)
plot(res_parametric_si, legend = FALSE)
dev.off()


