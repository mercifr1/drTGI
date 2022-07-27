###############################################################################

# TGI : Simulation SLD profiles from  ODE with mrgsolve 

# Scenario 1: unbalanced design matrix  
# 5 dose-cohorts with different number of participants per cohort and 
# different number of observations for each individual


# 31-07-2022

# Maia Muresan

################################################################################


# source ODE model 
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")
# source function to built Scenario 1
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_1_function.R")

#library (plyr)
#library( mrgsolve)

#' **Design matrix for Scenario 1**
#' -----------------------------------------------------------------------------


#' replicate the design matrix
replicate_unbalance_data<-replicate(2,{
  
  #' built the scenario 1 design matrix with function scenario_1() from Scenario_1_function
  unbalance_data<-scenario_1(doses = 10*c(1 / 6, 1 / 2, 1, 2, 3), 
                        nweek = 9, 
                        week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
                        n_N=5,
                        lambda_N=4,
                        lambda_0=4)
  list(unbalance_data) #' save data as list
})

#replicate_unbalance_data[[1]]

#' **Simulation**
#' ----------------------------------------------------------------------------
# simulate the SLD( i.e DV ) time profiles values

#' initialize a null list
unbalance_all_data<-list()
for(i in 1:length(replicate_unbalance_data)){
  #' simulation with mrgsim_d()
  unbalance_all_data[[i]]<-mrgsim_d(model_ode,as.data.frame(replicate_unbalance_data[[i]]),
                                    carry.out="DOSE,ALPHA,COH",
                                    end=-1,nocb=F)%>%
    dplyr::mutate(DV=RESP+EPS_1,
                  REP=i,
                  ID2=paste0(REP, "-", ID))
  
}


unbalance_design_matrix<-unbalance_all_data

#' transform into data frame 
final_unbalance_design_matrix <- plyr::ldply (unbalance_design_matrix, data.frame)


#''  remove patients that have 20 % increase from nadir
#'
final_unbalance_nadir<-final_unbalance_design_matrix %>%
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


#'plot the SLD profiles over time for one of the replicate
p_unbalance<-ggplot(filter(final_unbalance_nadir,REP==2),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  geom_point(size=2)+
  theme_bw()+
  ylab("SLD")+
  labs(color="DOSE",title="Scenario 1:Unbalanced design matrix",
       subtitle="Simulated SLD profiles")+
  scale_color_manual(values=c("#E3DB71","#E7298A","#1B9E77","#A6761D",  "#7570B3"  ))
  
p_unbalance
