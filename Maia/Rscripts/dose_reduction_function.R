
dose_reduction <- function(initial_data, cohort, proportion, reduction) {
  data_1 <- initial_data %>%
    filter(coh == cohort) %>%
    group_by(ID) %>%
    mutate(TIME_unif_reduc = round(runif(1, 7 * 12 / 356, 7 * 24 / 356), 3)) %>%
    ungroup()
  
  
  id_proportion<-proportion * length(unique(data1$ID))
  nsel<-ifelse(id_proportion-floor(id_proportion)<0.5,floor(id_proportion),ceiling(id_proportion))
  
  
  mysel_reduction <- sample(data_1$ID, size = nsel, replace = F)
  
  data_1 %>%
    mutate(
      IDflag_reduc = ifelse(ID %in% mysel_reduction, 1, 0),
      DOSE_red = ifelse(IDflag_reduc == 1 & TIME >= TIME_unif_reduc, reduction, DOSE)
    ) 
  
}


data_bind_reduc <- rbind(
  dose_reduction(initial_data = full_data, cohort = 3, proportion =  0.10, reduction = 5),
  dose_reduction(initial_data = full_data, cohort = 4, proportion =  0.25, reduction = 10),
  dose_reduction(initial_data = full_data, cohort = 5, proportion =  0.50, reduction = 20)
)



