################################################################################

# Fit TGI model:  a mixed-effect model defined with one ODE 

# Test distributions and the range of the parameters from the simulated data set
# Here we look only on  the design matrix from Scenario 1: unbalanced design matrix
# with flat dose-cohort

# 31-07-2022

# Maia Muresan

################################################################################


#
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Simulation_SLD_complete.R")
library(bayestestR)# for  Highest Density Interval (HDI)

# <design_flat>  data set from Simulation_SLD_complete.R
data_for_fit<-design_flat%>%filter(REP==1)




# Baseline
ggplot(data_for_fit, aes(x = BASE)) +
  geom_density(fill = "orange")

ggplot(data_for_fit, aes(x = ETA_4)) +
  geom_density(fill = "orange")
#credible interval based on the Highest Density Interval (HDI). 
#It will give us the range containing the 89% most probable effect values.
#Note that we will use 89% CIs instead of 95% CIs (as in the frequentist framework), 
#as the 89% level gives more stable results
hdi(data_for_fit$BASE, ci = 0.89)
# KG
ggplot(data_for_fit, aes(x = KG)) +
  geom_density(fill = "orange")
ggplot(data_for_fit, aes(x = ETA_1)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$KG, ci = 0.89)

#KS0
ggplot(data_for_fit, aes(x = KS0)) +
  geom_density(fill = "orange")
ggplot(data_for_fit, aes(x = ETA_2)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$KS0, ci = 0.89)

#GAMMA
ggplot(data_for_fit, aes(x = GAMMA)) +
  geom_density(fill = "orange")
ggplot(data_for_fit, aes(x = ETA_3)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$GAMMA, ci = 0.89)


# ALPHA
ggplot(data_for_fit, aes(x = ALPHA)) +
  geom_density(fill = "orange")
ggplot(data_for_fit, aes(x = ETA_5)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$ALPHA, ci = 0.89)

#RESPONSE variable
ggplot(data_for_fit, aes(x = DV)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$DV, ci = 0.89)


#error
ggplot(data_for_fit, aes(x = EPS_1)) +
  geom_density(fill = "orange")
#
hdi(data_for_fit$EPS_1, ci = 0.89)


#RESPONSE variable without the error
ggplot(data_for_fit, aes(x = (DV-EPS_1))) +
  geom_density(fill = "orange")
hdi(data_for_fit$DV-data_for_fit$EPS_1, ci = 0.89)
