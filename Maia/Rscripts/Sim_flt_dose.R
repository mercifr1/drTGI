library(mrgsolve)
library(purrr)
library(dplyr)


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

###'**Design matrix constructions**#--------------------------------------------

#'
pdtime<-rep(seq(0,1,0.1),times=5)
doses0<-rep(rep(10,11),times=5)


#'**design matrix with flate dose**
data0 <-
  data_frame(ID=rep(1:5,each=11)) %>%
  mutate(DOSE=doses0,TIME=pdtime)


###'**Design matrix for the analytical form**

##'Design matrix with flate dose 
##'
#' design matrix with evid 1 = a bolus or infusion dose
# data0_evid1 <-
#   data_frame(ID=rep(1:5,each=11)) %>%
#   mutate(DOSE=doses0,TIME=pdtime,CMT=0,EVID=1)

#' #' design matrix with  evid 2 = other type event, with solver reset
#' data0_evid2 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses0,TIME=pdtime,CMT=0,EVID=2)
#' 
#' 
#' #' design matrix with cmt=0 and evid 2 for DOSe=5 and 0 otherwise
#' data0_evid20 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses0,TIME=pdtime,CMT=0,EVID=ifelse(DOSE==5,2,0))


#'** simulation with flat dose**

#'  1. no cmt, no evid
set.seed(123)
out.ode<-model_ode %>%
  data_set(data0) %>%
  mrgsim(carry.out="DOSE")
out.ode.df<-as.data.frame(out.ode)


set.seed(123)
out.pred<-model_pred %>%
  
  data_set(data0) %>%
  mrgsim(carry.out="DOSE")
out.pred.df<-as.data.frame(out.pred)

data0%>%select(ID,DOSE,TIME)%>%
  mutate(Y_ODE_flt=out.ode.df$RESP,
         Y_PRED_flt=out.pred.df$RESP)%>%
  as.data.frame()%>%head(n=15)


#     ID DOSE TIME Y_ODE_flt Y_PRED_flt
# 1   1   10  0.0  83.90358   83.90358
# 2   1   10  0.1  81.45628   81.45628
# 3   1   10  0.2  79.63292   79.63292
# 4   1   10  0.3  78.35387   78.35387
# 5   1   10  0.4  77.55687   77.55687
# 6   1   10  0.5  77.19337   77.19337
# 7   1   10  0.6  77.22565   77.22565
# 8   1   10  0.7  77.62477   77.62476
# 9   1   10  0.8  78.36887   78.36887
# 10  1   10  0.9  79.44200   79.44200
# 11  1   10  1.0  80.83311   80.83311



###'

#' ##'Design matrix for flt dose
#' 

#' #' design matrix with cmt=1 and evid 1
#' data0_cmt1_evid1 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses0,TIME=pdtime,CMT=1,EVID=1)
#' 
#' 
#' #' design matrix with cmt=1 and evid 2
#' data0_cmt1_evid2 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses0,TIME=pdtime,CMT=1,EVID=2)
#' 
#' #' design matrix with cmt=1 and evid 2 for DOSE=5 and 0 otherwise
#' data0_cmt1_evid20 <-
#'   data_frame(ID=rep(1:5,each=11)) %>%
#'   mutate(DOSE=doses0,TIME=pdtime,CMT=1,EVID=ifelse(DOSE==5,2,1))

















