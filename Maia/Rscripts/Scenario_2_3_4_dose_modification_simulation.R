source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/functions.R")





unbalance_dose_modif<-replicate(10,{
  
#step 1: generate the initial data set  
 
  full_data<-scenario_1(doses = 10*c(1 / 6, 1 / 2, 1, 2, 3), 
                        nweek = 9,week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
                         n_N=5,lambda_N=4,lambda_0=4)

#step 2: based on the initial data set apply the dose reduction on the desired cohort

  data_bind_reduc <- rbind(
  dose_reduction(initial_data = full_data, cohort = 3, proportion =  0.10, reduction = 5,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356),
  dose_reduction(initial_data = full_data,cohort = 4, proportion =  0.25, reduction = 10,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356),
  dose_reduction(initial_data = full_data,cohort = 5, proportion =  0.50, reduction = 20,
                 n_unif_omis=1,min_unif_omis=7*12/356, max_unif_omis= 7*24/356)
)

#setp 3: based on the initial dataset we make the dose omission for the desired cohort

  data_bind_omis <- rbind(
  dose_omission(initial_data = full_data,cohort = 3,proportion = 0.10,omission = 0,
    n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
  ),
  dose_omission(initial_data = full_data,cohort = 4,proportion = 0.25,omission = 0,
                n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
  ),
  dose_omission(initial_data = full_data,cohort = 5,proportion = 0.50,omission = 0,
    n_unif_omis = 1, min_unif_omis = 7 * 12 / 356, max_unif_omis = 7 * 24 / 356
  )
)


#step 4: merge the dataset from the dose reduction and dose omission


  omis_reduc_merge<-merge(data_bind_reduc,data_bind_omis,
         by=c("ID","DOSE","COH","N_TIME","ID_OBS","TIME","EVID","CMT"),sort = F)



# step 5: compute the full dataset with all the cohorts and dose modifications

  
  sim_data_test<-full_data%>%
  filter(COH<=2)%>%
  mutate(TIME_unif_reduc=0,
         IDflag_reduc=0,
         DOSE_reduc=DOSE,
         TIME_unif_omis=0,
         IDflag_omis=0,
         TIMEflag_omis=0,
         DOSEflag_omis=0,
         DOSE_omis=DOSE)%>%
  rbind(omis_reduc_merge)%>%
  mutate(DOSEflag_omis=ifelse(DOSE_omis!=DOSE,1,0),
         DOSE_1=ifelse(IDflag_reduc==1 &DOSEflag_omis==0,DOSE_reduc,
                       ifelse(IDflag_reduc==0 & DOSEflag_omis==1,DOSE_omis,
                              ifelse(IDflag_reduc==1 & DOSEflag_omis==1,DOSE_omis,DOSE))))%>%
  dplyr::rename(DOSE_INITIAL=DOSE,DOSE=DOSE_1)%>%
  select(ID,DOSE,TIME,EVID,CMT,COH,N_TIME,ID_OBS)
  

  names(sim_data_test)<-toupper(names(sim_data_test))
  new_data<-sim_data_test

list(new_data)
})


#Simulation


unbalance_dose_modif_all<-list()
 for(i in 1:length(unbalance_dose_modif)){
   
   unbalance_dose_modif_all[[i]]<-mrgsim_d(model_ode,as.data.frame(unbalance_dose_modif[[i]]),carry.out="DOSE,ALPHA,COH",end=-1,nocb=F)%>%
   dplyr::mutate(DV=RESP+EPS_1,REP=i,ID2=paste0(REP, "-", ID))
   
 }
design<-unbalance_dose_modif_all

library (plyr)
scenario_2 <- plyr::ldply (design, data.frame)

scenario_2_nadir<-scenario_2 %>%
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
p_scenario2<-ggplot(filter(scenario_2_nadir,REP==2),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_scenario2






