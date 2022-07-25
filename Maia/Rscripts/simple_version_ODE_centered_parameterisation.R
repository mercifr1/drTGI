source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2_3_4_dose_modification_simulation.R")
library(tidyr)
library(dplyr)
library(rstan)
library(cmdstanr)

#' prepare the data for stan

# import the simulated dataset
data_for_fit <- design_flat %>%
  filter(REP == 1) %>%
  dplyr::mutate(TIME = ifelse(TIME < 0.1, 0, TIME), N_TIME = N_TIME) %>%
  select(ID, TIME, N_TIME, DOSE, DV, COH, ID_OBS)

#' buil a function to filter the *initial time point 0 (TRUE/FALSE)
data_prep <- function(data, no_0 = T) {
  if (no_0) {
    data_res <- filter(data, TIME != 0)
  } else {
    data_res <- data
  }
  data_res
}


gk_data <- data_prep(data_for_fit, no_0 = T)
nrow(gk_data)

data_gk <- list(
  T = nrow(gk_data), # nb. observations
  J = length(unique(gk_data$ID)), # nr. individuals
  sld = structure(gk_data$DV, .Dim = c(96, 1)), # response vector

  ts = gk_data %>% select(TIME) %>% pull(), # time points
  dose = gk_data$DOSE, # dose vector
  Nobs = gk_data %>% group_by(ID) %>%
    summarise("maxx" = max(N_TIME)) %>%
    select(maxx) %>% pull() # nb. of observations per individual
)





#' cmdstan compile
#'
mod <- cmdstan_model("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/simple_version_ODE.stan")

# check if teh syntax is correct
mod$check_syntax(pedantic = TRUE)


# fit
fit <- mod$sample(
  data = data_gk,
  iter_sampling = 3000,
  iter_warmup = 100, chains = 4,
  parallel_chains = 4,
  refresh = 10
)

# Results after sampling:

# All 4 chains finished successfully.
# Mean chain execution time: 100.7 seconds.
# Total execution time: 106.6 seconds.
#
# Warning: 11 of 12000 (0.0%) transitions ended with a divergence.
# See https://mc-stan.org/misc/warnings for details.
#
# Warning: 9741 of 12000 (81.0%) transitions hit the maximum treedepth limit of 10.
# See https://mc-stan.org/misc/warnings for details.
#
# Warning: 4 of 4 chains had an E-BFMI less than 0.2.
# See https://mc-stan.org/misc/warnings for details.
res <- fit$summary()



# plot
library(stringr) # useful for  Simple, Consistent Wrappers for Common String Operations
baseline_sld_data <- res %>% filter(str_detect(variable, "base")) # data with baseline values
sld_data <- res %>% filter(str_detect(variable, "DV")) # data with the actual sld values
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


fit$summary(variables = c("theta", "base", "lp__"), "mean", "sd")
fit$summary("theta", pr_lt_half = ~ mean(. <= 0.5))


# posterior draws
# default is a 3-D draws_array object from the posterior package
# iterations x chains x variables
draws_arr <- fit$draws() # or format="array"
str(draws_arr)
# draws x variables data frame
draws_df <- fit$draws(format = "df")
str(draws_df)
print(draws_df)
library(bayesplot)

mcmc_hist(fit$draws("theta"))


# this should be identical to draws_df created via draws(format = "df")
draws_df_2 <- as_draws_df(draws_arr)
identical(draws_df, draws_df_2)

#


#
plot_data <- res %>%
  filter(str_detect(variable, "DV")) %>%
  mutate(ID = gk_data$ID, TIME = gk_data$TIME, COH = gk_data$COH, DOSE = gk_data$DOSE)



ggplot(plot_data, aes(TIME, mean, group = ID, color = factor(DOSE))) +
  geom_line() +
  theme_bw() +
  # scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")
