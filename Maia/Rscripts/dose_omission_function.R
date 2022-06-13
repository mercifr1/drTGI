library(tidyr)
library(extraDistr)


dose_omission<-function(data, cohort, prop, omission){
  
  data1 <- data %>% #' data is from Scenario_2A
    dplyr::filter(coh == cohort) %>%
    dplyr::group_by(ID) %>%
    mutate(TIME_unif_omis = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
    ungroup()
  
  
  #' sample the 10% of subjects
   nsel<-floor(prop * length(unique(data1$ID)))

  mysel10_omission <- data1 %>%
    dplyr::sample_n(ID, size = nsel, replace = F)
  
   data1 %>%
    mutate(
      IDflag_omis = ifelse(ID %in% mysel10_omission, 1, 0), #' flag==1 for the specific subjects
      DOSE_omission = ifelse(IDflag_omis == 1 & TIME >= TIME_unif_omis & TIME < (TIME_unif_omis + 0.1), omission, DOSE)
    ) %>% #' one dose omission for those patients and at the specific time
    as.data.frame()
  
  
  
}
data_bind_omis<- rbind(
  dose_omission(data = data, cohort = 3, prop = 0.10, omission = 0.001),
  dose_omission(data = data, cohort = 4, prop = 0.25, omission = 0.001),
  dose_omission(data = data, cohort = 5, prop = 0.50, omission = 0.001)
)
