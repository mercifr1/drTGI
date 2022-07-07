source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2_3_4_dose_modification_simulation.R")

library(rstan)
library(cmdstanr)

dim(design_flat)
#' prepare the data for stan 
data_for_fit<-design_flat%>%
  filter(REP==1)%>%
  dplyr::mutate(TIME=ifelse(TIME<0.1,0,TIME))
         T
data_for_fit$timeInd<-seq(1:nrow(data_for_fit))

data <- list( N = nrow(data_for_fit),
              J = length(unique(data_for_fit$ID)),
              y = data_for_fit$DV,
              t0 = 0.0,
              timeInd= data_for_fit%>%select(timeInd)%>%pull(),
              ts=data_for_fit%>%select(TIME)%>%filter(TIME>0)%>%pull(),
              dose=data_for_fit$DOSE,
              T=data_for_fit %>% group_by(ID) %>% summarise("maxx" = max(N_TIME)) %>% select(maxx) %>% pull(),
              nTime = data_for_fit %>% group_by(ID) %>% summarise("maxx" = max(N_TIME)) %>% select(maxx) %>% pull()
)


#' cmdstan compile
mod<-cmdstan_model("Maia/Rscripts/TGI_ODE_stan_test_3.stan")



mod$print()
mod$exe_file()

fit <- mod$sample(
  data = data, 
  seed = 42L, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 0 # print update every 500 iters
)



#' Remarks:
#'Delay compilation:
#' in case we w ant to delay compilation when creating the CmdStanModel object
mod_delay<-cmdstan_model("Maia/Rscripts/TGI_ODE_stan_test_3.stan",compile = FALSE)
mod_delay$exe_file()
mod_delay$compile()
mod_delay$exe_file()

#'Pedantic check:
#' in case we want to check if there are any invalid syntax and also it enabled
#'    the check will also warn you about other potential issues in your as:
#' 1. distribution usages issues
#' 2. Unused parameter
#' 3. Large or small constant in a distr
#' 4. Control flow depends on a parameter: branching control (if/else)
#' 5. Variable is used before assignment
#' 6. strict or nonsensical parameter bounds

mod_pedantic<-cmdstan_model("Maia/Rscripts/TGI_ODE_stan_test_3.stan",pedantic=TRUE)

#' in case we want to turn on pedantic model separately from compilation:
mod_pedantic$check_syntax(pedantic = TRUE) 
#'or 

mod_pedantic <- cmdstan_model("Maia/Rscripts/TGI_ODE_stan_test_3.stan", compile = FALSE)
mod_pedantic$check_syntax(pedantic = TRUE)

#' Stan model variables
variables<-mod$variables()

names(variables)
names(variables$data)
names(variables$parameters)
names(variables$transformed_parameters)
names(variables$generated_quantities)


variables$data$J
variables$data$sigma
variables$parameters$kg
variables$transformed_parameters$theta





