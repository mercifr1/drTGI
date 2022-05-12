#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)
library(brms)

#'--------------------------ODE form: Simulation--------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model
#' 

## ODE without sigma
ode<-'
$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE=70,DOSE=10

$CMT RESP

$MAIN

double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0*(DOSE/10)*exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE=TVBASE*exp(ETA(4));
RESP_0=BASE;

double KS=KS0*exp(-GAMA*TIME);
ALAG_RESP=0.1;
$OMEGA  0.005 0.03 0.003 0.01


$ODE


dxdt_RESP = KG *RESP-KS*RESP;



$CAPTURE 
 RESP_0 BASE KS0 KS KG GAMA 
'
#' Compile the model

set.seed(123)
model_ode<-mcode("ODEtest",ode)





#' Construction of the Time matrix
#'----------------------------------------------------
nInd<-5 #number of individual

# buid a function to define the time 
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


#' Design matrix: one level dose for all patients for all time points
#'----------------------------------------------------


#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-c(10, 10, 10, 10, 10)

flt<-indf %>%
  mutate(TIME=ty,
         ID=as.numeric(as.character(ID)),DOSE=10,EVID=1,CMT=1,MDV=1)%>%
  subset(select=-c(ty))
flt2<-flt%>%mutate(DOSE=ifelse(tnum==5,5,DOSE),EVID=ifelse(tnum==5,0,1),CMT=1,MDV=1)
#flt1<-flt%>%mutate(DOSE=ifelse(tnum==5,5,DOSE))
head(flt,n=15)
head(flt2,n=15)


#' Simulation
#' ---------------------------------------------------

#' **1. For the flat dose**
#' 
#' Simulation without SIGMA
set.seed(123)
snd.flt_ODE<-mrgsim(model_ode, flt,carry.out="DOSE",end=-1,hmax=0.01) %>% 
  as.data.frame()
head(snd.flt_ODE,n=20)

set.seed(123)
snd.flt_ODE1<-mrgsim(model_ode, flt2,carry.out="DOSE") %>% 
  as.data.frame()
head(snd.flt_ODE1,n=20)
