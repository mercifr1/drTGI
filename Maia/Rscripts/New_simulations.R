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
DOSE=10, T=0

$OMEGA  0.005 0.03 0.003 0.01
$SIGMA 0.05

$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
capture RESP_0=BASE;
double KS =  KS0 * log(DOSE) * exp(-GAMA * T);
capture RESP = RESP_0*(exp((KG*T)-( KS0 * log(DOSE)/GAMA) + KS/GAMA)) + EPS(1);

$CAPTURE   BASE KS0 KS KG GAMA T
"


#' Compile the model
#' ------------------------------------------------------
set.seed(123)
mod <- mcode("maia", code)
mod
see(mod)
#' Time matrix
#'----------------------------------------------------
#'
#'

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5

#' buid a function to define the time 
funtime<-function(end,delta){
  ty<-round(seq(0,end,delta),3)
  p<-length(ty)
  return(data.frame(tnum=0:(p-1),ty))
}
# funtime(24,0.1)


indf<-tibble(ID=1:nInd) %>%
  split(.$ID) %>%
  map_dfr(., ~funtime(1,1/9), .id="ID")
indf


#' Design matrix: one level dose for all patients for all time points
#'----------------------------------------------------
#'
#'

#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-c(10, 10, 10, 10, 10)

flt<-indf %>%
  mutate(T=ty, 
         ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf["tnum"])+1))%>%
  subset(select=-c(ty))

head(flt,n=10)


#' Design matrix: 2 dose levels at a certain time point for all patients
#' ---------------------------------------------------

stepup0<-flt%>%mutate(DOSE=ifelse(flt$tnum==5,5,10))
head(stepup0,n=10)


#' Design matrix: 5 dose levels at a certain time point for all patients
#' ---------------------------------------------------
#' 
stepup1<-flt%>%mutate(DOSE=ifelse(tnum==2 | tnum==4|tnum==6| tnum==8,5,10))
head(stepup1,n=10)





#' Simulation 
#' ---------------------------------------------------
#' 
#' 

#' **1. For the flat dose**
#' 
set.seed(123)

snd.flt<-mrgsim_d(mod, as.data.frame(flt), carry.out="KS,DOSE,RESP_0") %>% 
  as.data.frame()%>%
  subset(select=-c(time))
head(snd.flt,n=10)

mrgsim_d(mod, as.data.frame(flt), carry.out="KS,DOSE,RESP_0")%>%plot(RESP~T)



#' **2. For the 2 level dose**
#' 
set.seed(123)

snd.stepup0<-mrgsim_d(mod, as.data.frame(stepup0), carry.out="KS,DOSE") %>% 
  as.data.frame()%>%
  subset(select=-c(time))

head(snd.stepup0,n=10)

set.seed(123)
mrgsim_d(mod, as.data.frame(stepup0), carry.out="KS,DOSE")%>%plot(RESP~T)



#' **2. For the 5 dose levels**
#' 
set.seed(123)

snd.stepup1<-mrgsim_d(mod, as.data.frame(stepup1), carry.out="KS,DOSE") %>% 
  as.data.frame()%>%
  subset(select=-c(time))

head(snd.stepup1,n=10)

set.seed(123)
mrgsim_d(mod, as.data.frame(stepup1), carry.out="KS,DOSE")%>%plot(RESP~T)



#'Fit with brms 
#'-----------------------


#'Prepare the data for fit
#'
#'
#'Flat Dose:
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
new_data<-snd.flt%>%
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMA=log(GAMA),lBAS=log(BASE))#%>%subset(select=-c(KS0,KG,GAMA,BASE))
names(new_data)

#' Define the model structure and parameters < random effects>
analytic <- bf(RESP ~ exp(lBAS) *  exp( exp(lKG) * T - ( exp(lKS0) *log(DOSE)/ exp(lGAMA) ) *( 1 - exp(exp(lGAMA)*T)))  ,
               nl = TRUE, 
               lBAS ~ 1+ (1|ID),
               lKG ~ 1 + (1 | ID),
               lKS0 ~ 1 + (1 | ID),
               lGAMA ~ 1 + (1 | ID))
analytic


#'** Define the priors**
#'
#'
#'
# check which  priors
getpr1<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * T - ( exp(lKS0) *log(DOSE)/ exp(lGAMA) ) *( 1 - exp(exp(lGAMA)*T)))  ,
          data = new_data, family = gaussian())


getpr2<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * T - ( exp(lKS0) *log(DOSE)/ exp(lGAMA) ) *( 1 - exp(exp(lGAMA)*T)))  ,
          data = new_data, family = lognormal())

getpr3<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * T - ( exp(lKS0) *log(DOSE)/ exp(lGAMA) ) *( 1 - exp(exp(lGAMA)*T)))  ,
                  data = new_data, family = student())

# # first try: it is not converging:max(Rhat) is  > 3
# 
# prior1<-set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lBAS")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lKG")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lKS0")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lGAMA")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lBAS")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lKG")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lKS0")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lGAMA")
# prior1
# 
# m01<-brm(analytic, data=new_data, prior=prior1, family=gaussian(),
#          
#          iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
# plot(m01)



# second try:It converge: Rhat=1 

#' Define the model structure and parameters < random effects>
analytic2 <- bf(RESP ~ exp(lBAS) *  exp( (exp(lKG)/1000) * T - ( (exp(lKS0)/1000) *log(DOSE)/ (exp(lGAMA)/1000) ) *( 1 - exp((exp(lGAMA)/1000)*T)))  ,
               nl = TRUE, 
               lBAS ~ 1+ (1|ID),
               lKG ~ 1 + (1 | ID),
               lKS0 ~ 1 + (1 | ID),
               lGAMA ~ 1 + (1 | ID))
analytic2

prior2<-set_prior("normal(4,1)",class="b",nlpar="lBAS")+
  set_prior("normal(0,1)",class="b",nlpar="lKG")+
  set_prior("normal(0,1)",class="b",nlpar="lKS0")+
  set_prior("normal(0,1)",class="b",nlpar="lGAMA")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lBAS")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKG")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKS0")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lGAMA")
prior2

m02<-brm(analytic2, data=new_data, prior=prior2, family=lognormal(link = "identity"),

         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m02)
plot(m02)

pairs(m02)

#' 2 Dose levels:
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
#' 
 

new_data_2<-snd.stepup0%>%
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMA=log(GAMA),lBAS=log(BASE))#%>%subset(select=-c(KS0,KG,GAMA,BASE))
names(new_data_2)
head(new_data_2)

m03<-brm(analytic2, data=new_data_2, prior=prior2, family=lognormal(link = "identity"),
         
         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m03)
plot(m03)


#'--------------------------ODE form: Simulation--------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model
#' 
ODEcode<-'

$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70

$PARAM@covariates 
DOSE=10,T=0

$CMT RESP

$MAIN

double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
capture BASE = TVBASE*exp(ETA(4));
RESP_0 = BASE;

$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0.05

$ODE
double KS=KS0 *log(DOSE)*exp( -GAMA * T);
dxdt_RESP = KG *RESP -   KS* RESP;


$CAPTURE 
KS  RESP_0 BASE
'

#' Compile the model
#' ------------------------------------------------------
set.seed(123)
ODEmod<-mcode("a1",ODEcode)
ODEmod
see(ODEmod)
param(ODEmod)

#' Time matrix
#'----------------------------------------------------
#'
#'

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5

#' buid a function to define the time ( funtime(,) from above)
# funtime(24,0.1)
funtime_ODE<-function(end,delta){
  time<-round(seq(0,end,delta),3)
  p<-length(time)
  return(data.frame(tnum=0:(p-1),time))
}

indf1<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~funtime_ODE(1,1/9), .id="ID")
head(indf1,n=10)

#' Design matrix: Flat dose
#'----------------------------------------------------

#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**
#' and to each ID corresponds a weight
doses<-c(10, 10, 10, 10, 10)

flt_ODE<-indf1 %>%
  mutate(T=time, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf1["tnum"])+1))
head(flt_ODE,n=10)


#' Design matrix: 2 level dose at a certain time point for all patients
#' ---------------------------------------------------
stepup0_ODE<-flt_ODE%>%mutate(DOSE=ifelse(flt_ODE$tnum==5,5,10))
head(stepup0_ODE)



#' Design matrix: 5 dose levels at a certain time point for all patients
#' ---------------------------------------------------
#' 
stepup1_ODE<-flt_ODE%>%mutate(DOSE=ifelse(tnum==2 | tnum==4|tnum==6| tnum==8,5,10))
head(stepup1_ODE)




#' Simulation
#' ---------------------------------------------------
#' 
#' **1. For the flat dose**
set.seed(123)
snd.flt_ODE<-mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="KS,KS0,DOSE,RESP_0,T,BASE") %>% 
  as.data.frame()
head(snd.flt_ODE)
mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="DOSE,T") %>%  plot(RESP~T)


#' **2. For the 2 level dose**
#' 
set.seed(123)
snd.stepup0_ODE<-mrgsim_d(ODEmod, as.data.frame(stepup0_ODE), carry.out="KS,DOSE,RESP_0,T,BASE") %>% 
  as.data.frame()

head(snd.stepup0_ODE,n=10)

mrgsim_d(ODEmod, as.data.frame(stepup0_ODE), carry.out="KS,DOSE,RESP_0,T")%>%plot(RESP~T)


#' **3. For the 5  dose levels**
#' 
set.seed(123)
snd.stepup1_ODE<-mrgsim_d(ODEmod, as.data.frame(stepup1_ODE), carry.out="KS,DOSE,RESP_0,T,BASE") %>% 
  as.data.frame()

head(snd.stepup1_ODE,n=10)

mrgsim_d(ODEmod, as.data.frame(stepup1_ODE), carry.out="KS,DOSE,RESP_0,T")%>%plot(RESP~T)


