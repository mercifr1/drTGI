###############################################################################

# TGI : Simulation SLD profiles from  ODE with mrgsolve 

# Scenario 0 balanced design matrix  
# 5 dose-cohorts with same number of observations for each individual

# 31-07-2022

# Maia Muresan

################################################################################
# source ODE model 
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")

#libraries
library (plyr)


# Make replications of the design matrix 


replicate_balance_data<-replicate(100,{
#dose-cohorts
doses = 10 * c(1 / 6, 1 / 2, 1, 2, 3)
# number of observations
nTA<-9 
#' dose administration at every 6 weeks
skel_6w<-7*c(0,6,12,18,24,30,36,42,48,52) 


# generate the number of individuals per cohort
individual_per_coh <- extraDistr::rtpois(n = 5, lambda = 4) 

# total number of individual per trail
total_number_indiv <- sum(individual_per_coh) 

# build a function to define the time 
funtime<-function(nTAm,skel){
  tweek<-jitter(skel,amount = 7)
  TIME<-round(tweek/365,3)
  return(data.frame(N_TIME=0:nTA,TIME))
}

# matrix with the equal number of observations per individual
id_observation_data<-tibble(ID = 1:total_number_indiv)%>%
  split(.$ID) %>%
  map_dfr(., ~funtime(nTA=nTA,skel=skel_6w),.id="ID")

#matrix with the number of patients per cohort 
id_cohort_data <- data.frame(DOSE = round(doses, 2),
                             COHORT = 1:length(doses), 
                             ID.PER.COH = individual_per_coh)%>%
  tidyr::uncount(., ID.PER.COH) %>%
  mutate(ID = 1:total_number_indiv)

#' merge the 2 data set by ID and keep all the columns from the first data set 
balance_data <- merge(x=id_cohort_data, y=id_observation_data, by = "ID", all.x = T)%>%
  dplyr::mutate(EVID=2,CMT=1)

list(balance_data)
})


# simulate the SLD( DV ) profile  values
balance_all_data<-list()

for(i in 1:length(replicate_balance_data)){
  
  balance_all_data[[i]]<-mrgsim_d(model_ode,as.data.frame(replicate_balance_data[[i]]),carry.out="DOSE,ALPHA,COH",end=-1,nocb=F)%>%
    dplyr::mutate(DV=RESP+EPS_1,REP=i,ID2=paste0(REP, "-", ID))
  
}

balance_design_matrix<-balance_all_data


final_balance_design_matrix <- plyr::ldply (balance_design_matrix, data.frame)



final_balance_nadir<-final_balance_design_matrix %>%
  group_by(ID)%>%
  mutate(MIN=min(DV),# find the minimum response variable SLD
         Tmin=TIME[which.min(DV)],# find the time corresponding to the minimum value
         FLAG_TIME=ifelse(TIME<=Tmin,0,1),# flag==1 for TIME > Tmin 
         MIN20=MIN+20*MIN/100,# compute the 20%incrise from nadir
         FLAG_NADIR20=ifelse(DV<=MIN20,0,1),
         FLAG_DVMIN=ifelse(FLAG_TIME==1 & FLAG_NADIR20==1,1,0)# flag==1 for DV> MIN20
  )%>%
  ungroup%>%
  filter(FLAG_DVMIN==0)%>%
  as.data.frame()

  
  #'plot
p_balance<-ggplot(filter(final_balance_nadir,REP==2),aes(TIME,DV,group=ID,color=factor(DOSE)))+
    geom_line(size=1)+
    theme_bw()+
    ylab("SLD")+labs(title="Scenario 0: Balance deisgn matrix with flat dose")
p_balance
  