################################################################################

# Fit TGI model:   Solving an simple version  ODE  for only 3 subject


#  Estimating System Parameters and Initial State with rstan

#  test only on design_flat: flat dose_cohort Scenario_1


# 31-07-2022

# Maia Muresan

################################################################################



library(rstan)
library(cmdstanr)
fit3<-stan (file="C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/simple_version_ODE_reff_two_subjects.stan",
           iter=1000,chains = 4)

results_3subjects<-summary(fit3)

