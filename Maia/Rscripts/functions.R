###############################################################################

# TGI : File with all the functions creates to built all scenarios
#     : Useful to generate the entire simulation

# function ti built the time matrix
# Scenario_1: initial data set unbalanced design matrix 
# Scenario_2: dose reduction
# Scenario_3: dose omission

# 31-07-2022

# Maia Muresan

################################################################################



#' build a function to define the time

funtime <- function(ny, week) {
  tweek <- jitter(week, amount = 7) #' add small amount of noise
  tyear <- round(tweek / 365, 3) #' transform time in year
  return(data.frame(N_TIME = 0:ny, tyear)) #' return data frame
}


# function to create scenario 1
scenario_1 <- function(doses, nweek, week,n_N,lambda_N, lambda_0) {
  #' doses: dose levels: a vector with the dose-cohorts
  #' nweek: dose administration at every 6 week : vector with the time( parameter for the for the funtime())
  #' week : #' total number of time  observations
  #' n_N  : integer: parameter to generate truncated Poisson distr. for the number of individual/cohort (here n_N=5 dose-cohorts)
  #' lambda_N: integer: parameter to generate the truncated Poisson distr.  
  #' lambda_0: integer: parameter to generate  the zero-truncated Poisson distr for the number of observations per individual
  #' NOTE: n_0 :=  total_number_indiv:  all the individual in the trail
  
  
  individual_per_coh <- extraDistr::rtpois(n = n_N, lambda = lambda_N) # generate the number of individuals per cohort
  total_number_indiv <- sum(individual_per_coh) # total number of individual per trail
  
  nb_obs_per_indiv <- actuar::rztpois(n = total_number_indiv, lambda = lambda_0) # generate number of observations per individual
  
  #' construct the data with the number of observations per individual
  id_obs_data <- tibble(ID = 1:total_number_indiv) %>%
    split(.$ID) %>%
    purrr::map_dfr(., ~ funtime(ny = nweek, week = week), .id = "ID") %>%
    mutate(ID_OBS = rep(nb_obs_per_indiv, each = max(.$N_TIME) + 1)) %>%
    dplyr::filter(N_TIME <= ID_OBS)
  
  #' construct the data with the  number of individual per cohort design matrix
  #' **To each cohort corresponds a different dose e.g. 5 cohort = 5 doses**
  #' 
  id_cohort_data <- data.frame(DOSE = round(doses, 2), COH = 1:length(doses), ID.PER.COH = individual_per_coh)%>%
    tidyr::uncount(., ID.PER.COH) %>%
    mutate(ID = 1:total_number_indiv)
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



# function to create scenario 2 with dose reduction on the desired cohorts
dose_reduction <- function(initial_data, cohort, proportion, reduction, n_unif_omis, min_unif_omis,max_unif_omis) {
  
  
  data_1 <- initial_data %>%
    filter(COH == cohort) %>%
    group_by(ID) %>%
    mutate(TIME_unif_reduc = round(x=runif(n=n_unif_omis, #generate 1 random uniform number per ID
                                           min=min_unif_omis, # minimum is the 2nd time dose administration
                                           max=max_unif_omis), # maximum is the 4th time dose administration
                                   digits=3)) %>%
    ungroup()
  
  #compute the proportion 
  id_proportion<-proportion * length(unique(data_1$ID))# proportion number of ID
  nsel<-ifelse(id_proportion-floor(id_proportion)<0.5,# round the proportion to a lower bound or upper bound based on the decimal
               floor(id_proportion),
               ceiling(id_proportion))
  
  #sample the proportion <nsel>% of subjects
  mysel_reduction <- sample(data_1$ID, size = nsel, replace = F)
  
  # compute the dose reduction for specific pacient
  data_1 %>%
    mutate(
      IDflag_reduc = ifelse(ID %in% mysel_reduction, 1, 0),# flag for the sampled ID
      
      # dose reduction starting from the Time larger than TIME_unif_reduc for the specific patients( IDflag_omis==1 )
      ID_DOSEflag_reduc=ifelse(IDflag_reduc == 1 & TIME >= TIME_unif_reduc,1,0),
      DOSE_reduc = ifelse(ID_DOSEflag_reduc==1, reduction, DOSE)
    ) 
  
}


# function to create scenario 3 with dose omission on the desired cohorts

dose_omission <- function(initial_data, cohort, proportion, omission, n_unif_omis, min_unif_omis, max_unif_omis) {
  
  
  data_1 <- initial_data %>% #' initial_data is from Scenario_1
    dplyr::filter(COH == cohort) %>%
    dplyr::group_by(ID) %>%
    mutate(TIME_unif_omis = round(x = runif(
      n = n_unif_omis,
      min = min_unif_omis,
      max = max_unif_omis
    ), digits = 3)) %>%
    ungroup()
  
  
  #' sample the proportion of subjects per cohort to receive the dose_omission
  id_proportion <- proportion * length(unique(data_1$ID))
  nsel <- ifelse(id_proportion - floor(id_proportion) < 0.5, floor(id_proportion), ceiling(id_proportion))
  
  mysel_omission <- sample(data_1$ID, size = nsel, replace = F)
  
  data_omis <- data_1 %>%
    group_by(ID) %>%
    mutate(
      IDflag_omis = ifelse(ID %in% mysel_omission, 1, 0), # flag for the suject in question
      TIMEflag_omis = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis, 1, 0),
      ID_DOSEflag_omis = {
        TIMEflag_omis == 1
      } %>%
        {
          . * !duplicated(.)
        },
      DOSE_omis = ifelse(ID_DOSEflag_omis == 1, omission, DOSE)
    ) %>%
    ungroup()
  
  for (i in unique(data_omis$ID)) {
    try1 <- filter(data_omis, ID == i & COH == 5) %>%
      pull(DOSE_omis)
    if (any(try1 == omission)) {
      place <- which(try1 == omission)
      if (length(try1) > place) {
        try1[place + 1] <- omission
      }
    }
    data_omis[data_omis$ID == i & data_omis$COH == 5, "DOSE_omis"] <- try1
  }
  
  data_omis%>%mutate(ID_DOSEflag_omis=ifelse(DOSE_omis!=DOSE,1,0))
}

