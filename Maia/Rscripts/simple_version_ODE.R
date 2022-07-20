source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2_3_4_dose_modification_simulation.R")

library(rstan)
library(cmdstanr)

#' prepare the data for stan 

# import the simulated dataset
data_for_fit<-design_flat%>%
  filter(REP==1)%>%
  dplyr::mutate(TIME=ifelse(TIME<0.1,0,TIME),N_TIME=N_TIME)%>%
  select(ID,TIME,N_TIME,DOSE,DV,BASE,KG,KS0,GAMMA,ALPHA,EPS_1,ID_OBS)

#' buil a function to filter the *initial time point 0 (TRUE/FALSE)
data_prep <- function(data, no_0 =T){
  if(no_0){
    data_res <- filter(data, TIME !=0)
  } else {
    data_res <- data
  }
  data_res
}


gk_data <- data_prep(data_for_fit, no_0 = T)
nrow(gk_data)

data_gk <- list(T = nrow(gk_data), # nb. observations
                J = length(unique(gk_data$ID)), # nr. individuals
                sld = structure(gk_data$DV,.Dim=c(96,1)), # response vector
                
                ts=gk_data%>%select(TIME)%>%pull(), # time points
                dose=gk_data$DOSE, # dose vector
                Nobs=gk_data %>% group_by(ID) %>% 
                  summarise("maxx" = max(N_TIME)) %>% 
                  select(maxx) %>% pull()# nb. of observations per individual
                
)


# initials values 
inits_1 <- list(tvkg=0.4,
                tvks0=0.4, 
                x_r=0.0,
                x_i=0,
                tvgamma= .2,
                tvalpha= .2,
                tvbase=70,
                
                omega_kg =.5,
                omega_ks0 =.5,
                omega_gamma=.5,
                omega_alpha=.5,
                omega_base=.5,
                sigma=0.5,
                eta_kg = rep(.5, 23),
                eta_ks0= rep(.5, 23),
                eta_gamma=rep(.5, 23),
                eta_alpha= rep(.5, 23),
                eta_base= rep(.5, 23)
                
                
)

# repeate the initials values to be the same as the number of chains in the sample Stan code
inits <- list(inits_1, inits_1, inits_1, inits_1)


#' cmdstan compile
#' 
mod<-cmdstan_model("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/simple_version_ODE.stan")

# check if teh syntax is correct
mod$check_syntax(pedantic = TRUE)


# fit
fit <- mod$sample(data = data_gk , 
                  iter_sampling = 3000,
                  iter_warmup = 1000,chains = 4,
                  parallel_chains = 4, 
                  refresh = 10)

# All 4 chains finished successfully.
# Mean chain execution time: 128.5 seconds.
# Total execution time: 132.5 seconds.
# Warning: 11571 of 12000 (96.0%) transitions hit the maximum treedepth limit of 10.
# See https://mc-stan.org/misc/warnings for details.
# 
# Warning: 4 of 4 chains had an E-BFMI less than 0.2.
# See https://mc-stan.org/misc/warnings for details.
res <- fit$summary()

