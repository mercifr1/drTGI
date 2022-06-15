

dose_reduction <- function(initial_data, cohort, proportion, reduction, n_unif_omis, min_unif_omis,max_unif_omis) {
  
  
  data_1 <- initial_data %>%
    filter(coh == cohort) %>%
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
      DOSE_reduc = ifelse(IDflag_reduc == 1 & TIME >= TIME_unif_reduc, reduction, DOSE)
    ) 
  
}

full_data<-scenario_1(doses = 10 * c(1 / 6, 1 / 2, 1, 2, 3), 
                      nweek = 9, 
                      week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
                      n_N=5,
                      lambda_N=4,
                      lambda_0=4)

data_bind_reduc <- rbind(
  dose_reduction(initial_data = full_data, 
                 cohort = 3, 
                 proportion =  0.10, 
                 reduction = 5,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356),
  dose_reduction(initial_data = full_data,
                 cohort = 4, 
                 proportion =  0.25, 
                 reduction = 10,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356),
  dose_reduction(initial_data = full_data,
                 cohort = 5, 
                 proportion =  0.50, 
                 reduction = 20,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356)
)


head(data_bind_reduc)

