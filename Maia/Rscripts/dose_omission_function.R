library(tidyr)
library(extraDistr)


dose_omission<-function(initial_data, cohort, proportion, omission){
  
  data_1 <- initial_data %>% #' data is from Scenario_2A
    dplyr::filter(coh == cohort) %>% 
    dplyr::group_by(ID) %>% 
    mutate(TIME_unif_omis = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>% 
    ungroup()
  
  
  #' sample the proportion of subjects per cohort to receive the dose_omision
  id_proportion<-proportion * length(unique(data_1$ID))
  nsel<-ifelse(id_proportion-floor(id_proportion)<0.5,floor(id_proportion),ceiling(id_proportion))
  
  mysel_omission <- data_1 %>%
    dplyr::sample_n(ID, size = nsel, replace = F)
  
  data_1 %>%
    mutate(
      IDflag_omis = ifelse(ID %in% mysel_omission, 1, 0),#flag for the suject in question
      TIMEflag_omis=ifelse( TIME>=TIME_unif_omis,1,0))%>% # flag for the time condition
    
    group_by(ID)%>%#  manipulation by ID group
    mutate(IDTIMEflag_omis = ifelse(IDflag_omis==1 & TIMEflag_omis==1,1, 0), #flag ID and TIME condition 
           ONE_REP=as.numeric(!duplicated(IDTIMEflag_omis)), # in case there are multiple conditions, select only the first equality
           DOSEflag_omis=ifelse(ONE_REP==1 & IDTIMEflag_omis==1,2,0),
           DOSE_omis=ifelse(DOSEflag_omis==2,omission,DOSE)#' one dose omission for those patients and at the specific time
    )%>%
    ungroup()
  
  
  
  
  
}

data_bind_omis<- rbind(
  dose_omission(initial_data = full_data, cohort = 3, proportion = 0.10, omission = 0.001),
  dose_omission(initial_data = full_data, cohort = 4, proportion = 0.25, omission = 0.001),
  dose_omission(initial_data = full_data, cohort = 5, proportion = 0.50, omission = 0.001)
)


