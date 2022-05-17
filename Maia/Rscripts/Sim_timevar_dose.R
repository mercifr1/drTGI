library(mrgsolve)
library(purrr)
library(dplyr)

#models

#'Analytical form without sigma
pred<- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70, DOSE=10

$OMEGA  0.005 0.03 0.003 0.01
$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
double RESP_0=BASE;


double KS =  KS0* exp(-GAMMA * TIME);
double RESP = RESP_0*exp((KG*TIME)-( KS0 *log(DOSE)/GAMMA)*(1-exp(-GAMMA*TIME)));

$CAPTURE    RESP RESP_0 BASE KS KS0 KG GAMMA EPS(1)
"
#' Compile the model
set.seed(123)
model_pred<-mcode("PRED",pred)


## ODE without sigma
ode<-'
$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70, DOSE=10

$CMT RESP


$MAIN
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE = TVBASE*exp(ETA(4));
RESP_0 = BASE;

$OMEGA  0.005 0.03 0.003 0.01

$ODE
double KS=KS0*exp( -GAMA * SOLVERTIME);
dxdt_RESP = KG *RESP -   KS*log(DOSE)*RESP;

$CAPTURE 
 RESP_0 BASE KS KS0 KG GAMA EPS(1) DOSE
'
#' Compile the model

set.seed(123)
model_ode<-mcode("ODE",ode)



#EVID is an event id indicator. mrgsolve recognized the following event IDs:

#0 = an observation record
#1 = a bolus or infusion dose
#2 = other type event, with solver reset
#3 = system reset
#4 = system reset and dose
#8 = replace


#'
pdtime<-rep(seq(0,1,0.1),times=5)
doses1<-rep(c(10,10,5,10,10,5,10,10,10,10,10),times=5)



#'**design matrix with time-varying dose**

data1 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses1,TIME=pdtime)

##'Design matrix with time-varying dose 
##'
#' design matrix for the analytical form with evid 1
data1_evid1 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses1,TIME=pdtime,CMT=0,EVID=1)

#' design matrix for the analytical form with  evid 2
data1_evid2 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses1,TIME=pdtime,CMT=0,EVID=2)


#' #' design matrix with cmt=0 and evid 2 for DOSe=5 and 0 otherwise
#' data1_evid20 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses1,TIME=pdtime,CMT=0,EVID=ifelse(DOSE==5,2,1))

##'Design matrix for time-varying dose
##'
#' design matrix for the ODE  with cmt=1 and evid 1
data1_cmt1_evid1 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses1,TIME=pdtime,CMT=1,EVID=1)


#' design matrix with cmt=1 and evid 2
data1_cmt1_evid2 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses1,TIME=pdtime,CMT=1,EVID=2)



#'** simulation with time-varying dose**

#'  1. no cmt, no evid
set.seed(123)
out.ode.tv<-model_ode %>%
  data_set(data1) %>%
  mrgsim(carry.out="DOSE")
out.ode.tv.df<-as.data.frame(out.ode_tv)


set.seed(123)
out.pred.tv<-model_pred %>%
  
  data_set(data1) %>%
  mrgsim(carry.out="DOSE")
out.pred.tv.df<-as.data.frame(out.pred.tv)

data1%>%select(ID,DOSE,TIME)%>%
  mutate(Y_ODE_tv=out.ode.tv.df$RESP,
         Y_PRED_tv=out.pred.tv.df$RESP)%>%
  as.data.frame()%>%head(n=15)
#    ID DOSE TIME Y_ODE_tv Y_PRED_tv
# 1   1   10  0.0 83.90358  83.90358
# 2   1   10  0.1 81.45628  81.45628
# 3   1    5  0.2 81.74429  84.08769
# 4   1   10  0.3 80.43132  78.35387
# 5   1   10  0.4 79.61319  77.55687
# 6   1    5  0.5 80.90288  87.19722
# 7   1   10  0.6 80.93672  77.22565
# 8   1   10  0.7 81.35501  77.62476
# 9   1   10  0.8 82.13487  78.36887
# 10  1   10  0.9 83.25957  79.44200
# 11  1   10  1.0 84.71753  80.83311
# 12  2   10  0.0 75.38454  75.38454
# 13  2   10  0.1 72.82724  72.82724
# 14  2    5  0.2 73.03295  75.36550
# 15  2   10  0.3 71.70652  69.66381


#'  1. cmt,  evid 1
set.seed(123)
out.ode.tv.evid1<-model_ode %>%
  data_set(data1_cmt1_evid1) %>%
  mrgsim(carry.out="DOSE,EVID,CMT",end=-1)
out.ode.tv.evid1.df<-as.data.frame(out.ode.tv.evid1)


set.seed(123)
out.pred.tv.evid1<-model_pred %>%
  
  data_set(data1_evid1) %>%
  mrgsim(carry.out="DOSE,EVID,CMT",end=-1)
out.pred.tv.evid1.df<-as.data.frame(out.pred.tv.evid1)

data1%>%select(ID,DOSE,TIME)%>%
  mutate(Y_ODE_tv1=out.ode.tv.evid1.df$RESP,
         Y_PRED_tv1=out.pred.tv.evid1.df$RESP,
         EVID_ODE=out.ode.tv.evid1.df$EVID,
         EVID_PRED=out.pred.tv.evid1.df$EVID)%>%
  as.data.frame()%>%head(n=15)


#'  2. cmt,  evid 2
set.seed(123)
out.ode.tv.evid2<-model_ode %>%
  data_set(data1_cmt1_evid2) %>%
  mrgsim(carry.out="DOSE,EVID,CMT",end=-1)
out.ode.tv.evid2.df<-as.data.frame(out.ode.tv.evid2)


set.seed(123)
out.pred.tv.evid2<-model_pred %>%
  
  data_set(data1_evid2) %>%
  mrgsim(carry.out="DOSE,EVID,CMT",end=-1)
out.pred.tv.evid2.df<-as.data.frame(out.pred.tv.evid2)

data1%>%select(ID,DOSE,TIME)%>%
  mutate(Y_ODE_tv2=out.ode.tv.evid2.df$RESP,
         Y_PRED_tv2=out.pred.tv.evid2.df$RESP,
         EVID_ODE=out.ode.tv.evid2.df$EVID,
         EVID_PRED=out.pred.tv.evid2.df$EVID)%>%
  as.data.frame()%>%head(n=15)
