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

# load the data 
data_all         = (read.csv("DATA/CORONA_TR.csv") %>% tbl_df())
cases_data_all   = data_all[,2];
deaths_data_all  = data_all[,4];
icu_data_all     = data_all[,6];
recov_data_all   = data_all[,7];
test_data_all    = data_all[,8];
popsize_data_all = 84180493;
allDates     = as.Date(parse_datetime(as.character(data_all$Date),"%d.%m.%y"));
date_data    = as.Date("2020-03-10");
date_end     = as.Date(tail(allDates, n=1));
data_start   = ymd(date_data)
data_end     = ymd(date_end)

# CORRECTED FOR THE NUMBER OF TESTS
test_data_all   = as.numeric(unlist(test_data_all))
newTestData     = daily_cases_data;
newTestData[47:length(daily_cases_data)] = max(test_data_all)*daily_cases_data[-(1:46)]/test_data_all[-(1:45)]
daily_cases_data   = newTestData
plot(as.incidence(daily_cases_data[-1], dates = allDates_agg[-1]))

m_relax_in = 0 
load(paste0(path2save,"/RDATA/T_modelTR_mrelax_",m_relax_in,".RData"))
     
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
png(filename=paste0(path2save,"/FIGS/RE_estimate_CORR.png"),
    width     = 6.25,
    height    = 6.25,
    units     = "in",
    res       = 1200)
plot(res_parametric_si, legend = FALSE)
dev.off()


