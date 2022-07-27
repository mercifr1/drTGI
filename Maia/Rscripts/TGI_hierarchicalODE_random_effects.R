################################################################################

# Fit TGI model:  a mixed-effect model defined with one ODE 

# Solving an  hierarchical ODE for time arrays of different size

# 31-07-2022

# Maia Muresan

################################################################################

source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2_3_4_dose_modification_simulation.R")
# how to install rstan: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(rstan)
#how to install cmdstan:https://mc-stan.org/docs/2_28/cmdstan-guide/cmdstan-installation.html
library(cmdstanr)

#' #####  Prepare the data for Stan #####

#' buil a function to filter the *initial time point 0 (TRUE/FALSE)
data_prep <- function(data, no_0 =T){
  if(no_0){
    data_res <- filter(data, TIME !=0)
  } else {
    data_res <- data
  }
  data_res
}


# import the simulated dataset <design_flat> from Scenario_2_3_4_dose_modification_simulation.R

# here we used only the simulation from the flat dose-cohort i.e  from scenario 1 design matrix
sim_data<-design_flat%>%
  filter(REP==1)%>%
  dplyr::mutate(TIME=ifelse(TIME<0.1,0,TIME),N_TIME=N_TIME)%>%
  select(ID,TIME,N_TIME,DOSE,DV,BASE,KG,KS0,GAMMA,ALPHA,EPS_1,ID_OBS)




sim_data_no_0 <- data_prep(sim_data, no_0 = T)
nrow(sim_data_no_t0)

data_stan <- list(N = nrow(sim_data_no_0),            # nb. observations
                uniqueT=nrow(sim_data_no_0),          # unique time points
                J = length(unique(sim_data_no_0$ID)), # nb. individuals
                
                y = sim_data_no_0$DV,                 # response vector
                t0 = 0.00,                            # initial time point
                ts=sim_data_no_0%>%select(TIME)%>%pull(), # time points
                dose=sim_data_no_0$DOSE,              # dose vector
                nTime=sim_data_no_0 %>% group_by(ID) %>% # nb. of observations per individual
                  summarise("maxx" = max(N_TIME)) %>% 
                  select(maxx) %>% pull()
)
                
# initials values 
inits_1 <- list(tvkg = .8,
                tvks0 = .4, 
                tvgamma = .4,
                tvalpha = 2,
                tvbase = 70,
                
                omega_kg = .5,
                omega_ks0 = .5,
                omega_gamma = .5,
                omega_alpha = .5,
                omega_base = .5,
                sigma   = .5,
                eta_kg = rep(.5, 23),
                eta_ks0 = rep(.5, 23),
                eta_gamma =rep(.5, 23),
                eta_alpha = rep(.5, 23),
                eta_base = rep(.5, 23)
)

# repeate the initials values to be the same as the number of chains in the sample Stan code
inits <- list(inits_1, inits_1, inits_1, inits_1)


#' cmdstan compile
#' 
mod<-cmdstan_model("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/TGI_ODE_reff.stan")

# check if the syntax is correct
mod$check_syntax(pedantic = TRUE)

# fit
fit <- mod$sample(data = data_gk , 
                  iter_sampling = 100,
                  iter_warmup = 50,
                  chains = 1,
                  parallel_chains = 1, 
                  refresh = 10)

res <- fit$summary()


dr <- fit$draws("tvkg", format = "draws_matrix")
saveRDS(dr, file = "exampleHowToSave.rds")


# fit
fit_more_chains <- mod$sample(data = data_gk , 
                  iter_sampling =3000,
                  iter_warmup = 2000,
                  chains = 4,
                  parallel_chains = 4, 
                  refresh = 500)
