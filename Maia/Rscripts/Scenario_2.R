source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")

#'Scenario 2 : unbalanced design matrix

#' -   dose-cohorts with different number of individual/cohort:
#'       - 5 cohorts
#'       - number of individual per cohort follows a truncated Poisson distr. (lambda_N=4)
#'       - obtain between 3 and 10 individuals per cohort.

#' -  each cohort have a different flat/constant dose=> 5 dose levels 
#' -  number of observations per individual is different and it follows a zero-truncated Poisson distr. (lambda_0=6)
#' -  lambda_0 reflect the short follow-up or early treatment discontinuation

library(purrr)

nInd<-25 #number of individual
skel_2<-7*c(0,6,12,18,24,30,36,42,48,52,58,64,70) # dose administration at every 6 week
nTA<-12 # number of observations

# build a function to define the time 
funtime<-function(nTAm,skel){
  tweek<-jitter(skel,amount = 7)
  ty<-round(tweek/365,3)
  return(data.frame(TAnum=0:nTA,ty))
}

indf.S2<-tibble(ID=1:nInd) %>%
  split(.$ID) %>%
  map_dfr(., ~funtime(nTA=nTA,skel=skel_2), .id="ID")
head(indf.S2,n=15)

#' Design matrix: different number of  individual per each cohort with Constant dose over time
#'----------------------------------------------------


#' **To each cohort corresponds a different dose e.g. 5 cohort = 5 doses**

doses<-10*c(1/6, 1/2, 1, 2, 3) 

#' Design matrix Scenario 2
flt.S2<-indf.S2 %>%
  mutate(TIME=ty,
         ID=as.numeric(as.character(ID)),
         DOSE=rep(round(doses,2), each=(max(indf.S2["TAnum"])+1)*5),
         TVKS0=ifelse(DOSE==1.67,0.8,ifelse(DOSE==30,0.3,0.4)),
         TVKG=ifelse(DOSE==1.67,0.3,ifelse(DOSE==30,0.8,0.6)),
         EVID=2,
         CMT=1)%>%
  subset(select=-c(ty))

head(flt.S2,n=15)


#' generate number of individual per cohort
#' install.packages("extraDistr")
library(extraDistr)# truncated-Poisson
set.seed(123)
ind.coh<-rtpois(5,4)
ind.coh


#' generate number of observations per individual 
#'install.packages("actuar")
library(actuar)#zero-truncated Poisson
set.seed(123)
ind.obs<-rztpois(25,6)
ind.obs



Cohorts<-flt.S2%>%
  mutate(coh=ifelse(ID<=3,1,
                    ifelse(ID<=9,2,
                           ifelse(ID<=12,3,
                                  ifelse(ID<=18,4,5)))))
Cohorts<-Cohorts%>%
  mutate(DOSE=ifelse(coh==1,1.67,
                     ifelse(coh==2,5,
                            ifelse(coh==3,10,
                                   ifelse(coh==4,20,30)))),
         OBS=rep(ind.obs,each=13))%>%
  filter(TAnum<OBS)


head(Cohorts,n=20)



#Simulation
set.seed(123)


ode.coh<-mrgsim_d(model_ode,as.data.frame(Cohorts),carry.out="DOSE,coh",end=-1,nocb=F)
ode.coh<-ode.coh%>%mutate(DV=RESP+EPS_1)
head(ode.coh,n=20)

p2<-ggplot(ode.coh, aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=0.8)+
  theme_bw()+
  xlab("time")+
  ylab("SLD")+
  guides(fill=guide_legend(title="DOSE"))
p2                  


