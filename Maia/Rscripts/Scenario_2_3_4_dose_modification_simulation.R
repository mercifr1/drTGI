source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/functions.R")


unbalance_dose_modif<-replicate(n=2,{
  
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

#step 3: based on the initial dataset we make the dose omission for the desired cohort

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



# step 5: compute the full dataset with all the cohorts and dose scenarios 1,2, 3 and 4

  
  sim_data_test<-full_data%>%
  filter(COH<=2)%>%
  mutate(TIME_unif_reduc=0,
         IDflag_reduc=0,
         ID_DOSEflag_reduc=0,
         DOSE_reduc=DOSE,
         TIME_unif_omis=0,
         IDflag_omis=0,
         TIMEflag_omis=0,
         ID_DOSEflag_omis=0,
         DOSE_omis=DOSE)%>%
  rbind(omis_reduc_merge)%>%
  mutate(ID_DOSEflag_reduc_omis=case_when(ID_DOSEflag_reduc==1 & ID_DOSEflag_omis==0~2,
                                          ID_DOSEflag_reduc==0 & ID_DOSEflag_omis==1~1,
                                          ID_DOSEflag_reduc==1 & ID_DOSEflag_omis==1~1,
                                 TRUE~0),
                  DOSE_reduc_omis=case_when(ID_DOSEflag_reduc_omis==2~DOSE_reduc,
                                            ID_DOSEflag_reduc_omis==1~DOSE_omis,
                                                                      TRUE~DOSE))
  
    

  names(sim_data_test)<-toupper(names(sim_data_test))
  new_data<-sim_data_test
  

list(new_data)

})




#Simulation


flat_dose_data<-list()
dose_redu_data<-list()
dose_omis_data<-list()
dose_reduc_omis_data<-list()
sim_flat_dose<-list()
sim_reduc_dose<-list()
sim_omis_dose<-list()
sim_reduc_omis_dose<-list()

#simulation of the flat dose

 for(i in 1:length(unbalance_dose_modif)){
   #simulation for the flat dose design matrix
   flat_dose_data[[i]]<-unbalance_dose_modif[[i]]%>%mutate(REP=i,
                                                          ID2=paste0(REP, "-", ID))
   sim_flat_dose[[i]]<-mrgsim_d(model_ode,as.data.frame(flat_dose_data[[i]]),carry.out="DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_REDUC,ID_DOSEFLAG_OMIS,ID_DOSEFLAG_REDUC_OMIS",end=-1,nocb=F)%>%
   dplyr::mutate(DV=RESP+EPS_1)
 }  

design_flat<-plyr::ldply (sim_flat_dose, data.frame)


design_flat_nadir<-design_flat %>%
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
p_design_flat<-ggplot(filter(design_flat_nadir,REP==1),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_design_flat




#simulation of the dose reduction
for(i in 1:length(unbalance_dose_modif)){
   #simulation for the  dose reduction design matrix
   dose_redu_data[[i]]<-unbalance_dose_modif[[i]]%>%mutate(DOSE=DOSE_REDUC,REP=i,
                                                           ID2=paste0(REP, "-", ID))
   sim_reduc_dose[[i]]<-mrgsim_d(model_ode,as.data.frame(dose_redu_data[[i]]),carry.out="DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_REDUC",end=-1,nocb=F)%>%
     dplyr::mutate(DV=RESP+EPS_1)
}  
design_reduction<-plyr::ldply(sim_reduc_dose,data.frame)

design_reduction_nadir<-design_reduction %>%
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
p_design_reduction<-ggplot(filter(design_reduction_nadir,REP==1),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_design_reduction

#simulation of the dose omission
for(i in 1:length(unbalance_dose_modif)){
   #simulation for the  dose omission design matrix
   dose_omis_data[[i]]<-unbalance_dose_modif[[i]]%>%mutate(DOSE=DOSE_OMIS,REP=i,
                                                           ID2=paste0(REP, "-", ID))
   sim_omis_dose[[i]]<-mrgsim_d(model_ode,as.data.frame(dose_omis_data[[i]]),carry.out="DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_OMIS",end=-1,nocb=F)%>%
     dplyr::mutate(DV=RESP+EPS_1)
} 
design_omission<-plyr::ldply(sim_omis_dose,data.frame)

design_omission_nadir<-design_omission %>%
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
p_design_omission<-ggplot(filter(design_omission_nadir,REP==1),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_design_omission



#simulation of the dose reduction and omission

for(i in 1:length(unbalance_dose_modif)){
   #simulation for the  dose omission design matrix
   dose_reduc_omis_data[[i]]<-unbalance_dose_modif[[i]]%>%mutate(DOSE=DOSE_REDUC_OMIS,REP=i,
                                                           ID2=paste0(REP, "-", ID))
   sim_reduc_omis_dose[[i]]<-mrgsim_d(model_ode,as.data.frame(dose_reduc_omis_data[[i]]),carry.out="DOSE,ALPHA,COH,REP,ID2,ID_DOSEFLAG_REDUC_OMIS",end=-1,nocb=F)%>%
     dplyr::mutate(DV=RESP+EPS_1)
   
}


design_reduction_omission<-plyr::ldply(sim_reduc_omis_dose,data.frame)

design_reduction_omission_nadir<-design_reduction_omission %>%
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
p_design_reduction_omission<-ggplot(filter(design_reduction_omission_nadir,REP==1),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_design_reduction_omission





