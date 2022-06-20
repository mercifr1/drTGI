#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr) #' useful for data manipulation
library(brms)

#'--------------------------ODE form: Simulation--------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model
#' 

## with ODE
ode<-'
$PARAM  TVKG =0.6, TVKS0=0.2, TVGAMMA=0.8, TVBASE =70, DOSE=10, TVALPHA=2

$OMEGA  0.005 0.03 0.003 0.01 0.04
$SIGMA  0.05
$CMT RESP
$MAIN
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE = TVBASE*exp(ETA(4));
double ALPHA = TVALPHA*exp(ETA(5));
RESP_0 = BASE;

$ODE
double KS=KS0 *exp( -GAMMA * SOLVERTIME);

if(DOSE<1) {dxdt_RESP = KG *RESP ;}
else{

dxdt_RESP = KG *RESP -   KS*ALPHA*log(DOSE)*RESP;}
$CAPTURE 
  RESP_0 ALPHA EPS(1) 
'
#' Compile the model

set.seed(123)
model_ode<-mcode("ODE",ode)


#' Construction of the Time matrix
#'----------------------------------------------------
nInd<-5 #number of individual

# build a function to define the time 
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


#' Design matrix: Constant dose over time
#'----------------------------------------------------


#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-c(10, 10, 10, 10, 10)

flto<-indf %>%
  mutate(TIME=ty,
         ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf["tnum"])+1),
         EVID=2,
         CMT=1)%>%
  subset(select=-c(ty))

head(flto,n=15)
#' Design matrix: Time-varying Dose
#'----------------------------------------------------

tvo<-flto%>%mutate(DOSE=ifelse(tnum==4,5,DOSE))



#' Simulation
#' ---------------------------------------------------

#' **1.  flat dose**
#' 
set.seed(123)
ode.flt<-mrgsim_d(model_ode, as.data.frame(flto),carry.out="DOSE",end=-1) %>% 
  as.data.frame()

ode.flt<-ode.flt%>%mutate(DV=RESP+EPS_1)
head(ode.flt,n=20)


ggplot(ode.flt, aes(TIME,DV,group=ID,color=factor(ID)))+
  geom_line()+
  theme_bw()+
  scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")

# set.seed(123)
# mrgsim_d(model_ode, as.data.frame(flto),carry.out="DOSE",end=-1)%>%plot(RESP~TIME)


#' **1. time-varying dose**
#' 
set.seed(123)
ode.tv<-mrgsim_d(model_ode, as.data.frame(tvo),carry.out="DOSE",end=-1,nocb=F) %>% 
  as.data.frame()
ode.tv<-ode.tv%>%mutate(DV=RESP+EPS_1)
head(ode.tv,n=20)


ggplot(ode.tv, aes(TIME,DV,group=ID,color=factor(ID)))+
  geom_line()+
  theme_bw()+
  scale_y_continuous(breaks = seq(50,90,5))+
  ylab("SLD")


# set.seed(123)
# mrgsim_d(model_ode, as.data.frame(tvo),end=-1)%>%plot(RESP~TIME)


