source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/functions.R")
#' Useful Libraries
#' -----------------------------------------------------------------------------
# install.packages("mrgsolve")
# library(mrgsolve) # for simulation

# install.packages("plyr")
# library(plyr )    # to transform list into dataframe



#' **Design matrix with scenarios 2, 3 and 4.**
#' ------------------------------------------------------------------------------
unbalance_dose_modif <- replicate(n = 10, { #' replicate the design matrix

  #' step 1: generate the initial data set with function <scenario_1> from functions.R
  #'

  full_data <- scenario_1(
    doses = 10 * c(1 / 6, 1 / 2, 1, 2, 3),
    nweek = 9, week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
    n_N = 5, lambda_N = 4, lambda_0 = 4
  )

  #' step 2: based on the initial data set apply the dose reduction on the desired cohort
  #' row bind the cohorts with the dose reduction
  #' use the function <dose_reduction> from functions.R

  data_bind_reduc <- rbind(
    dose_reduction(
      initial_data = full_data, cohort = 3, proportion = 0.10, reduction = 5,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    ),
    dose_reduction(
      initial_data = full_data, cohort = 4, proportion = 0.25, reduction = 10,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    ),
    dose_reduction(
      initial_data = full_data, cohort = 5, proportion = 0.50, reduction = 20,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    )
  )

  #' step 3: based on the initial data set apply the dose omission for the desired cohort
  #' row bind the cohorts with the dose omission
  #' use function <dose_omission> from functions.R
  #

  data_bind_omis <- rbind(
    dose_omission(
      initial_data = full_data, cohort = 3, proportion = 0.10, omission = 0,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    ),
    dose_omission(
      initial_data = full_data, cohort = 4, proportion = 0.25, omission = 0,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    ),
    dose_omission(
      initial_data = full_data, cohort = 5, proportion = 0.50, omission = 0,
      n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
    )
  )


  #' step 4: merge the resulting row bind data sets from the step 2: dose reduction and step 3: dose omission


  omis_reduc_merge <- merge(data_bind_reduc, data_bind_omis,
    by = c("ID", "DOSE", "COH", "N_TIME", "ID_OBS", "TIME", "EVID", "CMT"), sort = F
  )



  #' step 5: compute the full data set with all the cohorts and dose scenarios 1,2, 3 and 4
  #' row bind the full_data with dose reduction data set and dose omission data set

  sim_data_test <- full_data %>% #' initial data set
    dplyr::filter(COH <= 2) %>% #' filter the cohort with no dose modification
    #' add the corresponding columns to rowbind the dataset
    dplyr::mutate(
      TIME_unif_reduc = 0, #' Time for dose reduction
      IDflag_reduc = 0, #' IF flag for dose reduction
      ID_DOSEflag_reduc = 0, #' ID and TIME flag for dose_reduction
      DOSE_reduc = DOSE,
      TIME_unif_omis = 0, #' Time for dose omission
      IDflag_omis = 0, #' ID flag for dose omission
      TIMEflag_omis = 0, #' Time flag for dose omission
      ID_DOSEflag_omis = 0, #' ID and Time flag for dose omission
      DOSE_omis = DOSE
    ) %>%
    rbind(omis_reduc_merge) %>% #' row bind with the dose reduction and omission data set
    #' create  ID and DOSE flag for Scenario 4: accumulates events dose omission and dose reduction
    dplyr::mutate(
      ID_DOSEflag_reduc_omis = case_when(
        ID_DOSEflag_reduc == 1 & ID_DOSEflag_omis == 0 ~ 2, #' flag 2 for dose reduction when no dose omission
        ID_DOSEflag_reduc == 0 & ID_DOSEflag_omis == 1 ~ 1, #' flag 1 for dose omission when no dose reduction
        ID_DOSEflag_reduc == 1 & ID_DOSEflag_omis == 1 ~ 1, #' flag 1 for dose omission when there is a dose reduction
        TRUE ~ 0
      ), # otherwise 0
      DOSE_reduc_omis = case_when(
        ID_DOSEflag_reduc_omis == 2 ~ DOSE_reduc, #' replace the corresponding dose reduction amount when flag is 2
        ID_DOSEflag_reduc_omis == 1 ~ DOSE_omis, #' replace the corresponding dose omission mount when flag is 1
        TRUE ~ DOSE #' otherwise the corresponding flat dose
      )
    )



  names(sim_data_test) <- toupper(names(sim_data_test)) #' upper case for the columns name
  new_data <- sim_data_test #' design matrix with all the flag and dose modifications


  list(new_data) #' data set as list
})



#' **Simulation**
#' -----------------------------------------------------------------------------

#' simulation of the SLD profiles for all scenarios 2, 3 and 4


#' initialize null list for the design matrices for different scenarios
flat_dose_data <- list()
dose_redu_data <- list()
dose_omis_data <- list()
dose_reduc_omis_data <- list()

#' initialize null list for the simulation data for the corresponding scenarios
sim_flat_dose <- list()
sim_reduc_dose <- list()
sim_omis_dose <- list()
sim_reduc_omis_dose <- list()


#' Scenario_1 simulation of the flat dose (no dose modification)
#' -----------------------------------------------------------------------------

for (i in 1:length(unbalance_dose_modif)) {
  #' flat dose design matrix
  flat_dose_data[[i]] <- unbalance_dose_modif[[i]] %>%
    mutate(
      REP = i, #' corresponding replicate of the design matrix
      ID2 = paste0(REP, "-", ID)
    ) #' unique ID for each replicate
  #' simulation for the flat dose design matrix  with mrgsim_d()
  sim_flat_dose[[i]] <- mrgsim_d(model_ode, as.data.frame(flat_dose_data[[i]]),
    carry.out = "DOSE,COH,ID_OBS,N_TIME,REP,ID2",
    end = -1, nocb = F
  ) %>%
    dplyr::mutate(DV = RESP + EPS_1)
}
#' transform the simulated data set from list to data frame
design_flat <- plyr::ldply(sim_flat_dose, data.frame)

#' remove patients that have more than 20% increase from nadir
#'
design_flat_nadir <- design_flat %>%
  group_by(ID) %>% # group by subject
  mutate(
    MIN = min(DV), # find the minimum SLD response variable
    Tmin = TIME[which.min(DV)], # find the time corresponding to the minimum value Tmin
    FLAG_TIME = ifelse(TIME <= Tmin, 0, 1), # flag==1 for TIME > Tmin
    MIN20 = MIN + 20 * MIN / 100, # compute the 20% incrise from nadir
    FLAG_NADIR20 = ifelse(DV <= MIN20, 0, 1), # flag==1 for SLD smaller than possible 20 % increase nadir
    FLAG_DVMIN = ifelse(FLAG_TIME == 1 & FLAG_NADIR20 == 1, 1, 0) # flag==1 for DV > MIN20
  ) %>%
  ungroup() %>%
  filter(FLAG_DVMIN == 0) %>% #' select only the subjects less than 20% increse from nadir
  as.data.frame() #' save as data frame


#' plot the SLD profiles for one of the replicates
#'
p_design_flat <- ggplot(filter(design_flat_nadir, REP == 2), aes(TIME, DV, group = ID, color = factor(DOSE))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("SLD") +
  labs(color="DOSE",
    title = "Scenario 1: Unbalanced design matrix",
    subtitle = " SLD simulated profiles with  flat dose"
  )+scale_color_manual(values=c("#E3DB71","#E7298A","#1B9E77","#A6761D",  "#7570B3"  ))
p_design_flat




#' Scenario 2 simulation of the dose reduction
#' -----------------------------------------------------------------------------
for (i in 1:length(unbalance_dose_modif)) {
  #' complete design matrix with dose reduction
  dose_redu_data[[i]] <- unbalance_dose_modif[[i]] %>%
    mutate(DOSE = DOSE_REDUC, REP = i, ID2 = paste0(REP, "-", ID))

  #' simulation for the dose reduction
  sim_reduc_dose[[i]] <- mrgsim_d(model_ode, as.data.frame(dose_redu_data[[i]]),
    carry.out = "DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_REDUC",
    end = -1, nocb = F
  ) %>%
    dplyr::mutate(DV = RESP + EPS_1)
}

#' save matrix as data frame
design_reduction <- plyr::ldply(sim_reduc_dose, data.frame)


#' remove patients that have more than 20% increase from nadir
#'
design_reduction_nadir <- design_reduction %>%
  group_by(ID) %>%
  mutate(
    MIN = min(DV), # find the minimum response variable SLD
    Tmin = TIME[which.min(DV)], # find the time corresponding to the minimum value
    FLAG_TIME = ifelse(TIME <= Tmin, 0, 1), # flag==1 for TIME > Tmin
    MIN20 = MIN + 20 * MIN / 100, # compute the 20%incrise from nadir
    FLAG_NADIR20 = ifelse(DV <= MIN20, 0, 1),
    FLAG_DVMIN = ifelse(FLAG_TIME == 1 & FLAG_NADIR20 == 1, 1, 0) # flag==1 for DV> MIN20
  ) %>%
  ungroup() %>%
  filter(FLAG_DVMIN == 0) %>%
  as.data.frame()



#' plot the SLD profiles for one of the replicates
#'
p_design_reduction <- ggplot(filter(design_reduction_nadir, REP == 2), aes(TIME, DV, group = ID, color = factor(DOSE))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("SLD") +
  labs(color="DOSE",
    title = "Scenario 2: Unbalanced design matrix",
    subtitle = " SLD simulated profiles with  dose reduction"
  )+scale_color_manual(values=c("#E3DB71","#E7298A","#1B9E77","#A6761D",  "#7570B3"  ))
p_design_reduction

# "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"
# "#333BFF", "#CC6600", "#9633FF", "#E2FF33", "#E3DB71"
#' Scenario 3 simulation of the dose omission
#' -----------------------------------------------------------------------------
for (i in 1:length(unbalance_dose_modif)) {
  #' complete deisgn matrix  with dose omission
  dose_omis_data[[i]] <- unbalance_dose_modif[[i]] %>%
    mutate(
      DOSE = DOSE_OMIS, REP = i,
      ID2 = paste0(REP, "-", ID)
    )

  #' simulation for the  dose omission design matrix
  sim_omis_dose[[i]] <- mrgsim_d(model_ode, as.data.frame(dose_omis_data[[i]]),
    carry.out = "DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_OMIS",
    end = -1, nocb = F
  ) %>%
    dplyr::mutate(DV = RESP + EPS_1)
}

#' save simulated data set as data frame
design_omission <- plyr::ldply(sim_omis_dose, data.frame)


#' remove patients that have more than 20% increase from nadir
#'
design_omission_nadir <- design_omission %>%
  group_by(ID) %>%
  mutate(
    MIN = min(DV), # find the minimum response variable SLD
    Tmin = TIME[which.min(DV)], # find the time corresponding to the minimum value
    FLAG_TIME = ifelse(TIME <= Tmin, 0, 1), # flag==1 for TIME > Tmin
    MIN20 = MIN + 20 * MIN / 100, # compute the 20%incrise from nadir
    FLAG_NADIR20 = ifelse(DV <= MIN20, 0, 1),
    FLAG_DVMIN = ifelse(FLAG_TIME == 1 & FLAG_NADIR20 == 1, 1, 0) # flag==1 for DV> MIN20
  ) %>%
  ungroup() %>%
  filter(FLAG_DVMIN == 0) %>%
  as.data.frame()


#' plot the SLD profiles for one of the replicates
#'
p_design_omission <- ggplot(filter(design_omission_nadir, REP == 2), aes(TIME, DV, group = ID, color = factor(DOSE))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("SLD") +
  labs(color="DOSE",
    title = "Scenario 3: Unbalanced design matrix",
    subtitle = " SLD simulated profiles with dose omission"
  )+scale_color_manual(values=c("red","#E3DB71","#E7298A","#1B9E77","#A6761D",  "#7570B3"  ))
p_design_omission



#' Scenario 4 simulation of the cumulative events : Scenario 2 and Scenario 3
#' -----------------------------------------------------------------------------
for (i in 1:length(unbalance_dose_modif)) {
  #' complete design matrix with dose reduction and dose omission
  dose_reduc_omis_data[[i]] <- unbalance_dose_modif[[i]] %>%
    mutate(
      DOSE = DOSE_REDUC_OMIS, REP = i,
      ID2 = paste0(REP, "-", ID)
    )

  #' simulation
  sim_reduc_omis_dose[[i]] <- mrgsim_d(model_ode, as.data.frame(dose_reduc_omis_data[[i]]),
    carry.out = "DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_REDUC_OMIS",
    end = -1, nocb = F
  ) %>%
    dplyr::mutate(DV = RESP + EPS_1)
}

#' save simulated data set as data frame
design_reduction_omission <- plyr::ldply(sim_reduc_omis_dose, data.frame)


#' remove patients that have more than 20% increase from nadir
#'
design_reduction_omission_nadir <- design_reduction_omission %>%
  group_by(ID) %>%
  mutate(
    MIN = min(DV), # find the minimum response variable SLD
    Tmin = TIME[which.min(DV)], # find the time corresponding to the minimum value
    FLAG_TIME = ifelse(TIME <= Tmin, 0, 1), # flag==1 for TIME > Tmin
    MIN20 = MIN + 20 * MIN / 100, # compute the 20%incrise from nadir
    FLAG_NADIR20 = ifelse(DV <= MIN20, 0, 1),
    FLAG_DVMIN = ifelse(FLAG_TIME == 1 & FLAG_NADIR20 == 1, 1, 0) # flag==1 for DV> MIN20
  ) %>%
  ungroup() %>%
  filter(FLAG_DVMIN == 0) %>%
  as.data.frame()

#' plot the SLD profiles for one of the replicates
#'
p_design_reduction_omission <- ggplot(filter(design_reduction_omission_nadir, REP == 2), aes(TIME, DV, group = ID, color = factor(DOSE))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw() +
  ylab("SLD") +
  labs(color="DOSE",
    title = "Scenario 4: Unbalanced design matrix",
    subtitle = " SLD simulated profiles with  dose omission and dose reduction"
  )+scale_color_manual(values=c("red","#E3DB71","#E7298A","#1B9E77","#A6761D",  "#7570B3"  ))
p_design_reduction_omission



