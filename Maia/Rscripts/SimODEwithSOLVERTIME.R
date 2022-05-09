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
$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70

$PARAM@covariates 
DOSE=10

$CMT RESP

$MAIN
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE = TVBASE*exp(ETA(4));
RESP_0 = BASE;

$OMEGA  0.005 0.03 0.003 0.01


$ODE
double KS=KS0 *log(DOSE)*exp( -GAMA * SOLVERTIME);
dxdt_RESP = KG *RESP -   KS* RESP;


$CAPTURE 
 RESP_0 BASE KS KS0 KG GAMA EPS(1)
'
#' Compile the model

set.seed(123)
model_ode<-mcode("ODE",ode)


## ODE with SIGMA
ode_sigma<-'
$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMA=0.8, TVBASE =70

$PARAM@covariates 
DOSE=10

$CMT RESP

$MAIN
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));
double BASE = TVBASE*exp(ETA(4));
RESP_0 = BASE;

$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0.05

$ODE
double KS=KS0 *log(DOSE)*exp( -GAMA * SOLVERTIME);
dxdt_RESP = KG *RESP -   KS* RESP;

$TABLE
double DV=RESP+EPS(1);

$CAPTURE 
DV RESP_0 BASE KS KS0 KG GAMA EPS(1)
'
#' Compile the model

set.seed(123)
model_ode_sigma<-mcode("ODESIGMA",ode_sigma)





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
  mutate(time=ty, 
         ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=max(indf["tnum"])+1))%>%
  subset(select=-c(ty))

head(flt,n=15)

#' Simulation
#' ---------------------------------------------------

#' **1. For the flat dose**
#' 
#' Simulation without SIGMA
set.seed(123)
snd.flt_ODE<-mrgsim_d(model_ode, as.data.frame(flt)) %>% 
  as.data.frame()
head(snd.flt_ODE,n=20)

# Conclusion: <KS> is varying over time (NOTE: same values as the one in the Analytical simulation)
#           : <BASE> and <RESP_0> have the same values (NOTE: same values as the one in the Analytical simulation)
#           : <EPS(1)> are 0 overall
#           : <RESP> correspond to the correct value when the parameters values are replced in the formula

set.seed(123)
#plot the trend of RESP over time
mrgsim_d(model_ode, as.data.frame(flt)) %>%  plot(RESP~time)

set.seed(123)

#' **Simulation with SIGMA**
#' 
#' in the **<pred_sigma>** code we need to have some adjustement compared to **<pred>** code.
#' The **$ODE** is computing the solution without the error EPS(1), so to correspond to the mathematical/theoretical formula
#'    we shoul add to the solution of the dxdt_RESP i.e <RESP> the EPS(1), reason why we have :
#' $TABLE
#' double DV=REPS+EPS(1)
#' Thus, our solution and reponse variable is called : <DV> and not <RESP>

snd.flt_ODE_sigma<-mrgsim_d(model_ode_sigma, as.data.frame(flt)) %>% 
  as.data.frame()
head(snd.flt_ODE_sigma,n=20)

set.seed(123)
mrgsim_d(model_ode_sigma, as.data.frame(flt)) %>%  plot(DV~time)

# Conclusion: <KS> is varying over time (NOTE: same values as the one in the Analytical simulation)
#           : <BASE> and <RESP_0> have the same values (NOTE: same values as the one in the Analytical simulation, no <EPS(1)> in their calculation )
#           : <EPS(1)> are varying by individual over time
#           : <RESP> correspond to the correct value when the parameters values are replaced in the formula

