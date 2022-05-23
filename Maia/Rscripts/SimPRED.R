#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)
library(brms)


#'-------------------Analytically form: Simulation------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model

#'PRED
pred<- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70, DOSE=10

$OMEGA  0.005 0.03 0.003 0.01
$SIGMA 0.05

$PRED

double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
double RESP_0=BASE;


double KS =  KS0*exp(-GAMMA * TIME);

double RESP = RESP_0*exp((KG*TIME)-(KS0*log(DOSE)/GAMMA)*(1-exp(-GAMMA*TIME)));


$CAPTURE    RESP RESP_0 EPS(1)
"
#' Compile the model
set.seed(123)
model_pred<-mcode("PRED",pred)


#' Construction of Time matrix
#'----------------------------------------------------

nInd<-5 #number of individual

#' buid a function to define the time 
funtime<-function(end,delta){
  ty<-round(seq(0,end,delta),2)
  p<-length(ty)
  return(data.frame(tnum=0:(p-1),ty))
}
# funtime(24,0.1)

indf<-tibble(ID=1:nInd) %>%
  split(.$ID) %>%
  map_dfr(., ~funtime(1,1/10), .id="ID")
head(indf,n=15)


#' Design matrix: constant  dose over time
#'----------------------------------------------------

#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-c(10, 10, 10, 10, 10)

flt<-indf %>%
  mutate(TIME=ty, 
         ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf["tnum"])+1),
         EVID=2,
         CMT=0)%>%
  subset(select=-c(ty))

head(flt,n=15)

#' Design matrix: constant  dose over time
#'----------------------------------------------------

tv<-flt%>%mutate(DOSE=ifelse(tnum==4,5,10))
head(tv,n=15)

#' Simulation 
#' ---------------------------------------------------
#' 

#'Constant dose
set.seed(123)
snd.flt<-mrgsim(model_pred, flt, carry.out="DOSE",end=-1,output="df") 
snd.flt<-snd.flt%>%mutate(DV=RESP+EPS_1)
head(snd.flt)
#'Plot
#'
ggplot(snd.flt, aes(TIME,DV,group=ID,color=factor(ID)))+
  geom_line(size=1)+
  theme_bw()+
  scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")+
  title("Simulation Analytical form :Flat Dose")

#set.seed(123)
#mrgsim_d(model_pred, as.data.frame(flt),end=-1)%>%plot(RESP~TIME)







#'Time varying dose
#'
#' The PRED is not taking into account the time varying.
#' PRED is just a vector calculation: i.e take the values by each row
#' In our case, we want to capture the time-varying Dose
#' Thus: PRED should be computed recursively  
#' 
pred2 <- '
$PARAM TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70 ,DOSE=10
$OMEGA 0.005 0.03 0.003 0.01
$SIGMA 0.05

$PRED
double KG = TVKG * exp(ETA(1));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));

if(NEWIND <= 1)  {
  double RESP_0 = BASE;
  double TLAST = TIME;
  double KS0 = TVKS0*exp(ETA(2));
}

double DELT = TIME - TLAST;
double KS = KS0*exp(-GAMMA * DELT);
double RESP = RESP_0 * exp((KG*DELT)-(KS0*log(DOSE)/GAMMA)*(1-exp(-GAMMA*DELT)));


RESP_0 = RESP;
TLAST = TIME;
KS0 = KS;

$CAPTURE  RESP RESP_0 BASE EPS(1)
'

mod3 <- mcode("pred2", pred2, end = -1) 

set.seed(123)
snd.tv<-mrgsim_d(mod3, as.data.frame(tv),carry.out="DOSE") %>% 
  as.data.frame()

snd.tv<-snd.tv%>%mutate(DV=RESP+EPS_1)
head(snd.tv,n=20)

#'Plot
#'
ggplot(snd.tv, aes(TIME,DV,group=ID,color=factor(ID)))+
  geom_line(size=1)+
  theme_bw()+
  scale_y_continuous(breaks = seq(40,95,5))+
  ylab("SLD")

#' set.seed(123)
#' mrgsim_d(mod3, as.data.frame(tv),carry.out="DOSE") %>%plot(RESP~TIME)

