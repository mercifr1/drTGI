source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2A.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2_Dose_Reduction.R")


#' dose omission for coh 3-----------------------------------------------------
#'
df10_omission <- data %>% #' data is from Scenario_2A
  filter(coh == 3) %>%
  group_by(ID) %>%
  mutate(TIME_unif_omis = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
  ungroup()


#' sample the 10% of subjects
mysel10_omission <- df10_omission %>%
  dplyr::sample_n(ID, size = floor(0.10 * length(unique(ID))), replace = F)

df10_omission_new <- df10_omission %>%
  mutate(
    IDflag_omis = ifelse(ID %in% mysel10_omission, 1, 0), #' flag==1 for the specific subjects
    DOSE_omission = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis & TIME < (TIME_unif_omis + 0.1), 0.001, DOSE)
  ) %>% #' one dose omission for those patients and at the specific time
  as.data.frame()

#' length(unique(df10$ID))
#' round(ifelse((nt>=2&nt<=4),runif(8,7*12/356,7*24/356),TIME),3)


#' #' dose omission coh 4-------------------------------------------------------
#' 
df25_omission <- data %>% #' data is from Scenario_2A
  filter(coh == 4) %>%
  group_by(ID) %>%
  mutate(TIME_unif_omis = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
  ungroup()

#' sample the 25% of subjects
mysel25_omission <- df25_omission %>%
  dplyr::sample_n(ID, size = floor(0.25 * length(unique(ID))), replace = F)

df25_omission_new <- df25_omission %>%
  mutate(
    IDflag_omis = ifelse(ID %in% mysel25_omission, 1, 0), #' flag==1 for the specific subjects
    DOSE_omission = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis & TIME < (TIME_unif_omis + 0.1), 0.001, DOSE)
  ) %>% #' one dose omission  for those patients and at the specific time
  as.data.frame()


#' dose omission coh 5----------------------------------------------------------
#' 
df50_omission <- data %>% #' data is from Scenario_2A
  filter(coh == 5) %>%
  group_by(ID) %>%
  mutate(TIME_unif_omis = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
  ungroup()

#' sample the 50% of subjects
mysel50_omission <- df50_omission %>% dplyr::sample_n(ID, size = floor(0.50 * length(unique(ID))), replace = F)

df50_omission_new <- df50_omission %>%
  mutate(
    IDflag_omis = ifelse(ID %in% mysel50_omission, 1, 0), #' flag==1 for the specific subjects
    DOSE_omission = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis & TIME < (TIME_unif_omis + 0.2), 0.001, DOSE)
  ) %>% #' constraint :the dose omission always affects 2 dosing events in a row for the highest dose-cohort
  as.data.frame()


# nrow(df25[df25$coh == "4",])


#' Row Bind the data sets-------------------------------------------------------
#' 
#' row bind df10 and df25 and df_50
df10_25_omission <- rbind(df10_omission_new, df25_omission_new)
df10_25_50_omission <- rbind(df10_25_omission, df50_omission_new)

# select only the data with coh 1 and 2
data_omis <- data %>%
  filter(coh <= 2) %>%
  mutate(
    TIME_unif_omis = TIME,
    IDflag_omis = 0,
    DOSE_omission = DOSE
  )

#' row bind coh 1, 2,3 4, and 5
data_omission_final <- rbind(data_omis, df10_25_50_omission) %>% mutate(DOSE = DOSE_omission)

#' select the columns of interest
data_omission_final <- data_omission_final %>% select(ID, TIME, DOSE, nt, coh, EVID, CMT)

#' Simulation-------------------------------------------------------------------
#' 
set.seed(123)

ode.coh_omis <- mrgsim_d(model_ode, data_omission_final, carry.out = "DOSE,coh", end = -1, nocb = F)
ode.coh_omis <- ode.coh_omis %>%
  mutate(DV = RESP + EPS_1) %>%
  as.data.frame()


ode.coh_omis_nadir <- ode.coh_omis %>%
  group_by(ID) %>%
  mutate(
    MIN = min(DV),
    Tmin = TIME[which.min(DV)],
    FLAGt = ifelse(TIME <= Tmin, 0, 1),
    MIN20 = MIN + 20 * MIN / 100,
    FLAG20 = ifelse(DV <= MIN20, 0, 1)
  ) %>%
  ungroup() %>%
  as.data.frame()

#' plot-------------------------------------------------------------------------
p3 <- ggplot(ode.coh_omis, aes(TIME, DV, group = ID, color = factor(coh))) +
  geom_line() +
  theme_bw() +
  ylab("SLD")
p3
