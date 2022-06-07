source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")

#'Scenario 2_A : unbalanced design matrix. 
#' -  more general code for built the design matrix
#' -   dose-cohorts with different number of individual/cohort:
#'       - 5 cohorts
#'       - number of individual per cohort follows a truncated Poisson distr. (lambda_N=4)
#'       - obtain between 3 and 10 individuals per cohort.

#' -  each cohort have a different flat/constant dose=> 5 dose levels 
#' -  number of observations per individual is different and it follows a zero-truncated Poisson distr. (lambda_0=6)
#' -  lambda_0 reflect the short follow-up or early treatment discontinuation

library(purrr)# useful for data manipulation

nInd<-25 #number of individual
week<-7*c(0,6,12,18,24,30,36,42,48,52) # dose administration at every 6 week
nweek<-9# number of observations

# build a function to define the time 

funtime<-function(ny,week){
  tweek<-jitter(week,amount = 7)# add small amount of noise 
  tyear<-round(tweek/365,3) # transform time in year
  return(data.frame(nt=0:ny,tyear)) #return data frame
}


#construct the number of observations per individual
set.seed(123)
id_obs_data<-tibble(ID=1:nInd)%>%
  split(.$ID) %>%
  map_dfr(., ~funtime(ny=nweek,week=week), .id="ID")%>% 
  mutate(id_obs=rep(actuar::rztpois(n=nInd,lambda=4),each=10))%>%
  filter(nt<=id_obs)

head(id_obs_data,n=30)

#construct the number of individual per cohort design matrix
#' **To each cohort corresponds a different dose e.g. 5 cohort = 5 doses**

doses<-10*c(1/6, 1/2, 1, 2, 3)
set.seed(123)
id_cohort<-data.frame(DOSE=round(doses,2),coh=1:length(doses),ind.coh=extraDistr::rtpois(n=5,lambda=4))
id_cohort_data<-tidyr::uncount(id_cohort, ind.coh)%>%mutate(ID=1:25)
#uncount() : duplicating rows according to a weighting variable (or expression):performs the opposite to dplyr::count()
head(id_cohort_data)


#merge the 2 data set by ID and keep the DOSE  column from 
data<-merge(id_cohort_data,id_obs_data,by="ID",all.x=T)

data<-data%>%mutate(TIME=tyear,
                    ID=as.numeric(as.character(ID)),
                    TVKS0=ifelse(DOSE==1.67,0.8,ifelse(DOSE==30,0.3,0.4)),# change the parameters value 
                    TVKG=ifelse(DOSE==1.67,0.3,ifelse(DOSE==30,0.8,0.6)),# change the parameters value
                    EVID=2,
                    CMT=1)%>%
  subset(select=-c(tyear))






#Simulation
set.seed(123)

ode.coh<-mrgsim_d(model_ode,as.data.frame(data),carry.out="DOSE,coh",end=-1,nocb=F)
ode.coh<-ode.coh%>%mutate(DV=RESP+EPS_1)


ode.coh<-ode.coh%>%
  group_by(ID)%>%
  mutate(MIN=min(DV),
         Tmin=TIME[which.min(DV)],
         FLAGt=ifelse(TIME<=Tmin,0,1),
         MIN20=MIN+20*MIN/100,
         FLAG20=ifelse(DV<=MIN20,0,1)
  )%>%ungroup

p<-ggplot(ode.coh, aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line()+
  theme_bw()+
  scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")
p
