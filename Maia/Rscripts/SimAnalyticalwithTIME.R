#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)
library(brms)


#'-------------------Analytically form: Simulation------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model

code <- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70 
$PARAM@covariates
DOSE=10

$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0

$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
capture RESP_0=BASE;

double KS =  KS0 * log(DOSE) * exp(-GAMA * TIME);
capture RESP = RESP_0*(exp((KG*TIME)-( KS0 * log(DOSE)/GAMA) + KS/GAMA))+EPS(1);



$CAPTURE    RESP_0 KS KS0 KG GAMA EPS(1)
"

cod<- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70 
$PARAM@covariates
DOSE=10

$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0.05

$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
capture RESP_0=BASE;

double KS =  KS0 * log(DOSE) * exp(-GAMA * TIME);
capture RESP = RESP_0*(exp((KG*TIME)-( KS0 * log(DOSE)/GAMA) + KS/GAMA))+EPS(1);



$CAPTURE    RESP_0 KS KS0 KG GAMA EPS(1)
"


#' Compile the model
#' ------------------------------------------------------
set.seed(123)
mod <- mcode("maia", code)
mod
see(mod)
init(mod)
param(mod)
revar(mod)
omat(mod)
smat(mod)


set.seed(123)
model<-mcode("maia2",cod)
model
revar(model)
smat(model)
#' Time matrix
#'----------------------------------------------------
#'
#'

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5

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


#' Design matrix: one level dose for all patients for all time points
#'----------------------------------------------------
#'
#'

#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-c(10, 10, 10, 10, 10)

flt<-indf %>%
  mutate(time=ty, 
         ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf["tnum"])+1))%>%
  subset(select=-c(ty))
  
head(flt,n=15)

#' Design matrix: 2 dose levels at a certain time point for all patients
#' ---------------------------------------------------

stepup0<-flt%>%mutate(DOSE=ifelse(tnum==5,5,10.00))
head(stepup0,n=15)


#' Design matrix: 5 dose levels at a certain time point for all patients
#' ---------------------------------------------------
#' 
stepup1<-flt%>%mutate(DOSE=ifelse(tnum==2 | tnum==4|tnum==6| tnum==8,5,10))
head(stepup1,n=15)





#' Simulation 
#' ---------------------------------------------------
#' 
#' 

#' **1. For the flat dose with SIGMA=0**
#' 


set.seed(123)
#ssim data set
snd.flt<-mrgsim_d(mod, as.data.frame(flt), carry.out="KS,DOSE,RESP_0,EPS(1)") %>% 
  as.data.frame()
head(snd.flt,n=20)

set.seed(123)
# plot the data with SIGMA=0
mrgsim_d(mod, as.data.frame(flt), carry.out="KS,DOSE,RESP_0")%>%plot(RESP~time)


set.seed(123)
#plot the data when force the sigma to be 0 by using the function zero_re()
mrgsim_d(mod, as.data.frame(flt), carry.out="KS,DOSE,RESP_0",zero_re(sigma))%>%plot(RESP~time)

#' **1. For the flat dose wit SIGMA=0.05**
#' 
set.seed(123)

snd.flt_2<-mrgsim_d(model, as.data.frame(flt), carry.out="KS,DOSE,RESP_0,EPS(1)") %>% 
  as.data.frame()
head(snd.flt_2,n=20)
set.seed(123)

mrgsim_d(model, as.data.frame(flt), carry.out="KS,DOSE,RESP_0")%>%plot(RESP~time)

set.seed(123)

mrgsim_d(model, as.data.frame(flt), carry.out="KS,DOSE,RESP_0",zero_re(sigma))%>%plot(RESP~time)
set.seed(123)

model%>%zero_re(sigma)%>%data_set(flt)%>%mrgsim()%>%plot(RESP~time)
set.seed(123)

model%>%%>%data_set(flt)%>%mrgsim()%>%plot(RESP~time)


#' **2. For the 2 level dose**
#' 
set.seed(123)

snd.stepup0<-mrgsim_d(mod, as.data.frame(stepup0)) %>% 
  as.data.frame()

head(snd.stepup0,n=10)

set.seed(123)
mrgsim_d(mod, as.data.frame(stepup0))%>%plot(RESP~time)

# SIGMA=0.05
set.seed(123)

snd.stepup0_2<-mrgsim_d(model, as.data.frame(stepup0)) %>% 
  as.data.frame()

head(snd.stepup0_2,n=10)

set.seed(123)
mrgsim_d(model, as.data.frame(stepup0))%>%plot(RESP~time)



#' **3. For the 5 dose levels**
#' 
#' 

set.seed(123)

snd.stepup1<-mrgsim_d(mod, as.data.frame(stepup1), carry.out="KS,DOSE") %>% 
  as.data.frame()%>%
  subset(select=-c(time))

head(snd.stepup1,n=10)

set.seed(123)
mrgsim_d(mod, as.data.frame(stepup1), carry.out="KS,DOSE")%>%plot(RESP~T)


