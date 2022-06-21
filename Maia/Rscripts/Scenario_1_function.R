# Useful Libraries
#install.packages("tidyr")
library(tidyr)
library(dplyr)
#install.packages("purrr")
library(purrr)
#install.ackages("extraDistr")
library(extraDistr) # for extra distributions : i.e truncated Poisson
#install.packages("actuar")
library(actuar) # fro extra distributions: i.e zero-truncated Poisson


#' Scenario 1_function :  the unbalanced design matrix.
#'
#'
#'
#' -  built a function for the design matrix construction
#' 
#' -   dose-cohorts with different number of individual/cohort:
#'       - 5 dose-cohorts
#'       - number of individual per cohort follows a truncated Poisson distr. (lambda_N=4)
#'                 - obtain between 3 and 10 individuals per cohort.
#'       - total number of individual per trail = sum( number of indiv per cohort)
#'                                               =nb_id_coh1+nb_id_coh2+...+nb_id_coh5

#'       -  each cohort have a different flat/constant dose=> 5 dose levels
#'       
#' -  number of observations per individual is different and it follows a zero-truncated Poisson distr. (lambda_0=4)
#' 
#'       -  lambda_0 reflect the short follow-up or early treatment discontinuation
#'
#'


#' build a function to define the time

funtime <- function(ny, week) {
  tweek <- jitter(week, amount = 7) #' add small amount of noise
  tyear <- round(tweek / 365, 3) #' transform time in year
  return(data.frame(N_TIME = 0:ny, tyear)) #' return data frame
}

scenario_1 <- function(doses, nweek, week,n_N,lambda_N, lambda_0) {
  #' doses: dose levels: a vector with the dose-cohorts
  #' nweek: dose administration at every 6 week : vector with  time points ( parameter for the  funtime())
  #' week :  total number of time  observations: need to know length of vector nweek ( parameter for the  funtime())
  #' n_N  : integer: parameter to generate truncated Poisson distr. rtpois() for the number of individual/cohort (here n_N=5 dose-cohorts)
  #' lambda_N: integer: parameter to generate the truncated Poisson distr. rtpois() 
  #' lambda_0: integer: parameter to generate  the zero-truncated Poisson distr rztpois() for the number of observations per individual
  #' NOTE: n_0 :=  total_number_indiv:  all the individual in the trail
  
  #' generate random numbers 
  individual_per_coh <- extraDistr::rtpois(n = n_N, lambda = lambda_N) # generate the number of individuals per cohort
  total_number_indiv <- sum(individual_per_coh) # total number of individual per trail

  nb_obs_per_indiv <- actuar::rztpois(n = total_number_indiv, lambda = lambda_0) # generate number of observations per individual

  #' construct the data with the number of observations per individual
  id_obs_data <- tibble(ID = 1:total_number_indiv) %>%
    split(.$ID) %>%
    purrr::map_dfr(., ~ funtime(ny = nweek, week = week), .id = "ID") %>%
    dplyr::mutate(ID_OBS = rep(nb_obs_per_indiv, each = max(.$N_TIME) + 1)) %>%
    dplyr::filter(N_TIME <= ID_OBS)

  #' construct the data with the  number of individual per cohort design matrix
  #' **To each cohort corresponds a different dose e.g. 5 cohort = 5 doses**
  #' 
  id_cohort_data <- data.frame(DOSE = round(doses, 2), COH = 1:length(doses), ID.PER.COH = individual_per_coh)%>%
    tidyr::uncount(., ID.PER.COH) %>%
    dplyr::mutate(ID = 1:total_number_indiv)
  #' NOTE: uncount() : duplicating rows according to a weighting variable (or expression):performs the opposite to dplyr::count()

  #' merge the 2 data set by ID and keep all the columns from the first data set 
  data <- merge(x=id_cohort_data, y=id_obs_data, by = "ID", all.x = T)

  data %>%
    mutate(
      TIME = tyear,
      ID = as.numeric(as.character(ID)),
      EVID = 2,
      CMT = 1
    ) %>%
    subset(select = -c(tyear))
}



full_data<-scenario_1(doses = 10 * c(1 / 6, 1 / 2, 1, 2, 3), 
                      nweek = 9, 
                      week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
                      n_N=5,
                      lambda_N=4,
                      lambda_0=4)



