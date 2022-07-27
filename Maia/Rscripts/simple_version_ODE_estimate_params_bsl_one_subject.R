################################################################################

# Fit TGI model:   Solving an simple version  ODE  for only 1 subject


#  Estimating System Parameters and Initial State with rstan

# test only on design_flat: flat dose_cohort Scenrio_1


# 31-07-2022

# Maia Muresan

################################################################################

source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Simulation_SLD_complete.R")


library(rstan)

simulated_data<-design_flat%>%
  filter(REP==1)%>%
  dplyr::mutate(TIME=ifelse(TIME<0.1,0,TIME),N_TIME=N_TIME)%>%
  select(ID,TIME,N_TIME,ID_OBS,DOSE,DV,BASE, KG, KS0, GAMMA, ALPHA)

head(simulated_data$TIME,n=16)
head(simulated_data$DV,n=16)





T <- 20
t0 <- 0
theta <- c(0.6,0.15,0.8,1.6)
dose<-rep(1.67,times=T)
ts <- seq(0.1,1,length.out=20)
 # c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
y_hat <- structure(head(simulated_data$DV,n=20),.Dim=c(20,1))
  ##structure(c(0.617979, -0.365649, -0.816086, -0.513888, 0.0250277, 0.618301,
            #   0.307184, 0.0610513, -0.415887, -0.386049, 0.0257084, 0.368119, 0.445961,
            #   0.141085, -0.143833, -0.321638, -0.209449, 0.243388, 0.324905, -0.0382953,
            #   -0.773039, -0.608271, 0.0815231, 0.44686, 0.541102, 0.209472, -0.419926,
            #   -0.609051, -0.28225, 0.29107, 0.358185, 0.223183, -0.162348, -0.4618,
            #   -0.0397883, 0.11546, 0.446885, 0.252464, -0.0585741, -0.129213),
            # .Dim = c(20, 2))



data<-list(T=T,t0=t0,ts=ts,theta=theta, y_hat=y_hat,dose=dose)


fit1<-stan (file="C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/simple_version_ODE_ref_one_subject.stan",
           iter=1000,chains = 1)

results_1subject<-summary(fit1)# is a list 
summary_params<-results_1subject[[1]]
params<-results_1subject[[2]]

summaries_fit1<-summary(fit)$summary

# Save an object to a file
sum_data<-saveRDS(summaries_fit1, file = "simple_version_ODE_ref_one_subject_fit.rds")
# Restore the object
readRDS(file = "simple_version_ODE_ref_one_subject_fit.rds")
