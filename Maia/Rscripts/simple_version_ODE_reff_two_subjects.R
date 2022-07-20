library(rstan)

simulated_data<-design_flat%>%
  filter(REP==1)%>%
  dplyr::mutate(TIME=ifelse(TIME<0.1,0,TIME),N_TIME=N_TIME)%>%
  select(ID,TIME,N_TIME,ID_OBS,DOSE,DV,BASE, KG, KS0, GAMMA, ALPHA)

head(simulated_data$TIME,n=16)

head(simulated_data$DV,n=16)











fit_2<-stan (file="C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/simple_version_ODE_reff_two_subjects.stan",
           iter=1000,chains = 1)

results_2<-summary(fit)

summary(fit)
head(data_for_fit$DV,n=8)
head(design_flat$TIME,n=8)
head(design_flat$DOSE,n=8)