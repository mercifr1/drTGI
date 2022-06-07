source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_2A.R")


#'Scenario with dose-reduction in cohort 3, 4 and 5 

#'step1: make the changes on the cohorts in question independent on for the other
#' a) cohort 3: 10% subjects receive dose reduction from 10*1 mg/Q3W to 10*1/2 mg/Q3W
#' b) cohort 4: 25% subjects receive dose reduction from 10*2 mg/Q3W to 10*1 mg/Q3W mg/Q3W
#' c) cohort 5: 50% subjects receive dose reduction from 10*3 mg/Q3W to 10*2 mg/Q3W mg/Q3W


#'step2: (row) bind all the data sets


#'Step1:
#' **NOTE:** for the desired patients, after dose reduction, the dose amount remain the same (dose not come back to the initial amount) 
#' 
#' a) create the dose reduction for 10% of subjects on coh 3
#' 
set.seed(123)
df10 <- data %>%
  filter(coh == 3) %>% 
  group_by(ID) %>%
  mutate(TIME_unif = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
  ungroup() #' select coh 3 from data; create a uniform distr. variable for the time reduction

#' sample the 10% of subjects
mysel10 <- sample(df10$ID, size = ceiling(0.10 * length(unique(df10$ID))), replace = F)
#'
df10_new <- df10 %>%
  mutate(
    IDflag = ifelse(ID %in% mysel10, 1, 0), #' flag==1 for the specific subjects
    DOSE2 = ifelse(IDflag == 1 & TIME >= TIME_unif, 5, DOSE)
  ) %>% #' reduce the dose for those patients and at the specific time
  as.data.frame()



#' b) create the dose reduction for 25% of subjects on coh 4
set.seed(123)
df25 <- data %>%
  filter(coh == 4) %>% #' select coh 4
  mutate(TIME_unif = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) #' crate the uniform distr for the time reduction

#' sample the 25% of subjects
set.seed(123)
mysel25 <- unique(sample(df25$ID, size = ceiling(0.25 * length(unique(df25$ID))), replace = F))

#'
df25_new <- df25 %>%
  mutate(
    IDflag = ifelse(ID %in% mysel25, 1, 0), #' flag==1 for the specific subjects
    DOSE2 = ifelse(IDflag == 1 & TIME >= TIME_unif, 10, DOSE)
  ) %>% #' reduce the dose for those subjects at the specific time
  as.data.frame()


#' c) create the dose reduction for 50% of subjects on coh 5
set.seed(123)
df50 <- data %>%
  filter(coh == 5) %>% # select coh 5
  mutate(TIME_unif = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) #' create the uniform distr for the time reduction

#' sample the 50% of subjects
set.seed(123)
mysel50 <- unique(sample(df50$ID, size = ceiling(0.50 * length(unique(df50$ID))), replace = F))
#
df50_new <- df50 %>%
  mutate(
    IDflag = ifelse(ID %in% mysel50, 1, 0), #' flag==1 for the specific subjects
    DOSE2 = ifelse(IDflag == 1 & TIME >= TIME_unif, 20, DOSE)
  ) %>% #' reduce the dose for those subjects at the specific time
  as.data.frame()


# row bind df10 and df25
df10_25 <- rbind(df10_new, df25_new)
df10_25_50 <- rbind(df10_25, df50_new)



#' step2': filter the original data for coh 1 and 2
data_coh12 <- data %>%
  filter(coh <= 2) %>%
  mutate(
    TIME_unif = TIME,
    IDflag = 0,
    DOSE2 = DOSE
  )
data_final <- rbind(data_coh12, df10_25_50) %>% mutate(DOSE = DOSE2)
data_final <- data_final %>% select(ID, TIME, DOSE, nt, coh, EVID, CMT)

#' Simulation
set.seed(123)

ode.coh_reduction <- mrgsim_d(model_ode, data_final, carry.out = "DOSE,coh", end = -1, nocb = F)
ode.coh_reduction <- ode.coh_reduction %>%
  mutate(DV = RESP + EPS_1) %>%
  as.data.frame()


ode.dose_reduction_nadir<- ode.coh_reduction %>%
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

#' plot
p1 <- ggplot(ode.coh_reduction, aes(TIME, DV, group = ID, color = factor(coh))) +
  geom_line() +
  theme_bw() +
  ylab("SLD")
p1

