source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Scenario_1_function.R")
library (plyr)



replicate_unbalance_data<-replicate(1000,{
  
  unbalance_data<-scenario_1(doses = 10*c(1 / 6, 1 / 2, 1, 2, 3), 
                        nweek = 9, 
                        week = 7 * c(0, 6, 12, 18, 24, 30, 36, 42, 48, 52),
                        n_N=5,
                        lambda_N=4,
                        lambda_0=4)
  list(unbalance_data)
})

#replicate_unbalance_data[[1]]


# simulate the SLD( DV ) values
unbalance_all_data<-list()
for(i in 1:length(replicate_unbalance_data)){
  
  unbalance_all_data[[i]]<-mrgsim_d(model_ode,as.data.frame(replicate_unbalance_data[[i]]),carry.out="DOSE,ALPHA,COH",end=-1,nocb=F)%>%
    dplyr::mutate(DV=RESP+EPS_1,REP=i,ID2=paste0(REP, "-", ID))
  
}

unbalance_design_matrix<-unbalance_all_data


final_unbalance_design_matrix <- plyr::ldply (unbalance_design_matrix, data.frame)



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


#'plot
p_unbalance<-ggplot(filter(final_unbalance_nadir,REP==20),aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=1)+
  theme_bw()+
  ylab("SLD")
p_unbalance
