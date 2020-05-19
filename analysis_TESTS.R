# Setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("setup.R")
library(anytime)
library(readr)
library(splines)
library(forecast)


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
data_control_1 = ymd("2020-03-18") # first lockdown (assuming 2 lockdowns)
date_control_1 = ymd("2020-03-18")
data_control_2 = ymd("2020-04-03") # second lockdown (assuming 2 lockdowns)
date_control_2 = ymd("2020-04-03")
date_control_3 = ymd("2020-03-21") # mid lockdown (assuming 1 lockdown)
data_control_3 = ymd("2020-03-21")
date_relax     = ymd("2020-05-11") # relaxation (if any)

test_data_all   = as.numeric(unlist(test_data_all))
test_data_all_N = test_data_all/max(test_data_all)
recov_data_all  = as.numeric(unlist(recov_data_all))
recov_data_all_N= recov_data_all/max(recov_data_all)

x_all = seq(from=1,to=length(allDates))


# fit to lognormal function
fit_1  = nls(test_data_all_N ~ K*(1/(x_all*sig*sqrt(2*pi)))*exp(-((log(x_all)-mu)/(sqrt(2)*sig))^2), data = data.frame(x_all, test_data_all_N),
             start = c(K = 70, mu = 4, sig=0.5))
test_fit_vec_1=as.numeric(coef(fit_1))

png(filename=paste0(path2save,"/FIGS/TESTS_estimate.png"),
    width     = 6.25,
    height    = 5.25,
    units     = "in",
    res       = 1200)
plot(allDates,fitted(fit_1),type="l", col="blue", lwd=2, xlab="gun", ylab="normalize test sayisi",ylim=range(fitted(fit_1),test_data_all_N),panel.first = grid())
lines(allDates,test_data_all_N,pch=19,col="black",type="b", lty=2,panel.first = grid())
legend(allDates[2],0.8,legend=c("Lognormal fonksiyon", "Veri"), col=c("blue", "black"), lty=1:2, cex=0.8)
dev.off()


