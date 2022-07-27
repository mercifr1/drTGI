###############################################################################

# TGI : Built function  to create Scenario 2 design matrix

# Unbalanced design matrix with dose reduction
# Dose reduction  on the highest dose-cohort at a specific time point

# 31-07-2022

# Maia Muresan

################################################################################
library(tidyr)
library(extraDistr)

# Built a function to create the design matrix with dose reduction on highest dose-cohort


# We set to 10%, 25% and 50% the proportion of patients with dose omission in
# the 3rd, 4th and 5th dose-cohorts respectively, leaving the other cohorts unaffected

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
      DOSE_reduc = ifelse(IDflag_reduc == 1 & TIME >= TIME_unif_reduc, reduction, DOSE)
    ) 
  
}



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

