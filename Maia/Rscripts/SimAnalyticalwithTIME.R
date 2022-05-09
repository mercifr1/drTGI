#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)
library(brms)


#'-------------------Analytically form: Simulation------------------------------
#'------------------------------------------------------------------------------

#'Define the structure of the model

#'Without SIGMA
#NOTE : we can still capture EPS(1), because by default $SIGMA matrix in 0

pred<- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70 
$PARAM@covariates
DOSE=10

$OMEGA  0.005 0.03 0.003 0.01

$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
double RESP_0=BASE;

double KS =  KS0 * log(DOSE) * exp(-GAMMA * TIME);
double RESP = RESP_0*(exp((KG*TIME)-( KS0 * log(DOSE)/GAMMA) + KS/GAMMA));

$CAPTURE    RESP RESP_0 BASE KS KS0 KG GAMMA EPS(1)
"
#' Compile the model
set.seed(123)
model_pred<-mcode("PRED",pred)


#' With SIGMA
# NOTE: we add  +EPS(1) to the RESP, compared to the <pred> model

pred_sigma<- "
$PARAM    TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70 
$PARAM@covariates
DOSE=10

$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0.05

$PRED
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE=TVBASE * exp(ETA(4));
double RESP_0=BASE;

double KS =  KS0 * log(DOSE) * exp(-GAMMA * TIME);
double RESP = RESP_0*(exp((KG*TIME)-( KS0 * log(DOSE)/GAMMA) + KS/GAMMA))+EPS(1);

$CAPTURE    RESP RESP_0 BASE KS KS0 KG GAMMA EPS(1)
"

#' Compile the model
set.seed(123)
model_pred_sigma<-mcode("PREDSIGMA",pred_sigma)



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
#' 

#' Simulation for pred ( without SIGMA)
set.seed(123)
snd.flt<-mrgsim_d(model_pred, as.data.frame(flt), carry.out="DOSE") %>% 
  as.data.frame()
head(snd.flt,n=20)

# Conclusion: <KS> in varying over <time> ( what we expected)
#           : <BASE> and <RESP_0> have the same value ( what we expected)
#           : <EPS(1)> is 0 overall
#           : <RESP> values are correct when we replace all the parameters inside the formula
# Ex: for time=0.4: RESP_0*(exp((KG*TIME)-( KS0 * log(DOSE)/GAMA) + KS/GAMA))+EPS(1):

#     77.55688=83.90358*(exp((0.6429033 *0.4)-( 0.4236797  * log(10)/0.7705246) + 0.7168025/0.7705246))+0

# plot
set.seed(123)
mrgsim_d(model_pred, as.data.frame(flt), carry.out="DOSE") %>% plot(RESP~time)



#'Simulation for pred_sigma ( with SIGMA 0.05)
set.seed(123)
snd.flt_sigma<-mrgsim_d(model_pred_sigma, as.data.frame(flt),carry.out="DOSE") %>% 
  as.data.frame()
head(snd.flt_sigma,n=20)


#Conclusion: KS is the same as in the previous data ( correct)
#          : BASE and RESP_0 remains the same as before ( correct )
#          : EPS(1) is changing by individual over <time>
#          : <RESP> values are correct when we replace the parameters + <EPS(1)> 
#            inside the formula, and they are different from the <RESP> from the 
#            simulation without SIGMA
# Ex: for time=0.4: 
#     77.45759=83.90358 *(exp((0.6429033 *0.4)-( 0.4236797  * log(10)/0.7705246 ) + 0.7168025/0.7705246 ))+0.222903976
set.seed(123)
mrgsim_d(model_pred_sigma, as.data.frame(flt)) %>% plot(RESP~time)
