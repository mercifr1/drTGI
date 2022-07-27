################################################################################

# Fit TGI model:  a mixed-effect model defined with one ODE 

#  Stan code for solving  a simple version of hierarchical ODE
#                      with time array of different sizes 

# 31-07-2022

# Maia Muresan

################################################################################


#Remark: We only test the fit on the flat dose-cohort design matrix scenario 1 
# source 
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Simulation_SLD_complete.R")

#libraries
library(tidyr)
library(dplyr)
# How to install rstan : https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(rstan)
# How to install cmdstan: https://mc-stan.org/cmdstanr/
library(cmdstanr)
#install.packages("straingr")
library(stringr) # useful for  Simple, Consistent Wrappers for Common String Operations




#' Prepare the data for stan

#' built a function to filter the *initial time point 0 (TRUE/FALSE)
data_prep <- function(data, no_0 = T) {
  if (no_0) {
    data_res <- filter(data, TIME != 0)
  } else {
    data_res <- data
  }
  data_res
}


# import the first repetition simulated data set <design_flat> for flat dose-cohorts
data_for_fit <- design_flat %>%
  filter(REP == 1) %>%
  dplyr::mutate(TIME = ifelse(TIME < 0.1, 0, TIME), N_TIME = N_TIME) %>%
  select(ID, TIME, N_TIME, DOSE, DV, COH, ID_OBS)

data_wtho_t0<- data_prep(data_for_fit, no_0 = T)
#nrow(s_data)


# data for stan 
data_stan <- list(
  T = nrow(data_wtho_t0), # nb. observations
  J = length(unique(data_wtho_t0$ID)), # nr. individuals
  sld = structure(data_wtho_t0$DV, .Dim = c(96, 1)), # response vector
  
  ts = data_wtho_t0 %>% select(TIME) %>% pull(), # time points
  dose = data_wtho_t0$DOSE, # dose vector
  Nobs = data_wtho_t0 %>% group_by(ID) %>%
    summarise("maxx" = max(N_TIME)) %>%
    select(maxx) %>% pull() # nb. of observations per individual
)





#' cmdstan compile
#'
mod_stan <- cmdstan_model("simple_version_ODE_centered_parameterisation.stan")

# check if the syntax is correct
mod_stan$check_syntax(pedantic = TRUE)


# sampling
fit_simple_version <- mod_stan$sample(
  data = data_gk,
  iter_sampling = 3000,
  iter_warmup = 100, chains = 4,
  parallel_chains = 4,
  refresh = 10
)

# Results after sampling:

# All 4 chains finished successfully.
# Mean chain execution time: 113.6 seconds.
# Total execution time: 116.7 seconds.
# 
# Warning: 12 of 12000 (0.0%) transitions ended with a divergence.
# See https://mc-stan.org/misc/warnings for details.
# 
# Warning: 9489 of 12000 (79.0%) transitions hit the maximum treedepth limit of 10.
# See https://mc-stan.org/misc/warnings for details.
# 
# Warning: 4 of 4 chains had an E-BFMI less than 0.2.
# See https://mc-stan.org/misc/warnings for details.


# summary 
res_simple_version <- fit_simple_version$summary()



# plot sld 



baseline_sld_data <- res_simple_version %>% filter(str_detect(variable, "base")) # data with baseline values
sld_data <- res_simple_version %>% filter(str_detect(variable, "DV")) # data with the actual sld values
full_data <- rbind(baseline_sld_data, sld_data) # row bind baseline and sld dataset

ID_new <- c(1:nrow(baseline_sld_data), gk_data$ID) # define the ID column
TIME_new <- c(rep(0, times = 23), gk_data$TIME) # define TIME column
DOSE_new <- data_for_fit$DOSE # define the DOSE column

full_data <- full_data %>%
  mutate(ID = ID_new, TIME = TIME_new) %>%
  arrange(ID) %>%
  mutate(DOSE = DOSE_new, COH = data_for_fit$COH)

ggplot(full_data, aes(TIME, mean, group = ID, color = factor(COH))) +
  geom_line(size = 0.8) +
  geom_point(size = 0.8) +
  theme_bw() +
  # scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")

#


fit_simple_version$summary(variables = c("kg","ks0","gamma","alpha", "base", "lp__"), "mean", "sd")
fit_simple_version$summary("theta", pr_lt_half = ~ mean(. <= 0.5))


# posterior draws
# default is a 3-D draws_array object from the posterior package
# iterations x chains x variables
draws_arr <- fit_simple_version$draws() # or format="array"
str(draws_arr)
# draws x variables data frame
draws_df <- fit_simple_version$draws(format = "df")
str(draws_df)
print(draws_df)



