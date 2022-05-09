
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimAnalyticalwithTIME.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODEwithSOLVERTIME.R")
library(dplyr)
library(ggplot2)




ODE<-snd.flt_ODE_sigma
head(ODE,n=20)

ODE<-ODE%>%subset(select=-c(RESP))%>%rename(RESP=DV)
head(ODE)

ANALYTIC<-snd.flt_sigma



head(ANALYTIC,n=20)
head(ODE,n=20)



set.seed(123)

par(mfrow=c(2,1))
set.seed(123)
mrgsim_d(model_pred_sigma, as.data.frame(flt)) %>% plot(RESP~time)
set.seed(123)
mrgsim_d(model_ode_sigma, as.data.frame(flt)) %>% plot(DV~time)
