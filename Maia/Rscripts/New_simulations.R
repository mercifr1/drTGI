

#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)



library(brms)


#-------------------Analytically form: Simulation------------------------------
#------------------------------------------------------------------------------



#### First try-----------------------------------------------------------------

set.seed(123)


### Implemet the structure of the model

code <- "

$PARAM  TVKG = 0.001, TVKS0=0.2, TVGAMA=0.9, DOSE=10, T=0

$OMEGA 0.005 0.03 0.003

$SIGMA 0.1

$PRED
double Y0 = 70;

double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));

double KS =  KS0 * log(DOSE) * exp(-GAMA * T);

capture Y = Y0*(exp((KG*T)-( KS0 * log(DOSE)/GAMA)+KS))+EPS(1);


$CAPTURE  KS0 KS KG GAMA TIME 
"


### Compile the model code
set.seed(123)
mod <- mcode("maia", code)
mod


# illus_DOSE<-data.frame(DOSE=c(10,10 , 10, 5, 7,10,10,10,10,10),T=seq(0,1,1/9)) %>%
#   
#   mutate(KG = 0.001 * exp(0.005),
#          KS0 = 0.2 * exp(0.03),
#          GAMA = 0.9 * exp(0.003),
#          KS = KS0 * DOSE * exp(-GAMA * TIME),
#          Y = 70*(exp((KG*TIME)-(KS0 * DOSE/GAMA)+KS))+0.1
#                          )
# 
# g01_ybreaks<-c(10,10 , 10, 5, 7,10,10,10,10,10)
# g01_ylabs<-as.character(c(10,10 , 10, 5, 7,10,10,10,10,10))
# g01<-ggplot(illus_DOSE, aes(DOSE, Y))+
#   geom_point()+
#   geom_line()+
#   scale_y_continuous("Y (unit)", breaks=unique(illus_DOSE$DOSE))+
#   scale_x_continuous("DOSE (mg Q3W iv)", breaks=g01_ybreaks, labels=g01_ylabs)+
#   theme_minimal()+
#   theme(panel.grid.minor=element_blank())
# g01
# 





#' Time matrix
#'----------------------------------------------------
#'
#'
set.seed(123)



func_weeks<-function(nTAm, skel){
  rTAdys<-round(skel,2)
          
  TAweek<-round(rTAdys, 3)              #' Round to 3 digit years
  return(data.frame(TAnum=0:nTA, TAweek))
} 




#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5
nTA<-9
skel<-seq(0,1,1/9)

set.seed(123)
indf0<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~func_weeks(nTAm = nTA,skel = skel), .id="ID")

indf0


#' Design matrix: one level dose for all patients for all time points
#'----------------------------------------------------

#' To each ID corresponds a different dose e.g. 5 ID = 5 doses

set.seed(123)
doses<-c(10, 10, 10, 10, 10)

flt<-indf0 %>%
  mutate(T=TAweek, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=nTA+1))

flt



#' Simulation 
#' ---------------------------------------------------
set.seed(123)
snd.flt<-mrgsim_d(mod, as.data.frame(flt), carry.out="T,KS,KS0,DOSE") %>% as.data.frame()%>%subset(select = -c(time,TIME))
snd.flt

mrgsim_d(mod, as.data.frame(flt), carry.out="T") %>%  plot(Y+KS~T)


############### Fit with brms ##################################################


#prepare the data for fitting

new_data<-snd.flt%>%mutate(lKG=log(KG),lKS0=log(KS0),lGAMA=log(GAMA))%>%subset(select=-c(KS0,KG,GAMA))
names(new_data)
analytic <- bf(Y ~ exp(4.25) * (exp((exp(lKG) * T) - (exp(lKS0) *log(DOSE)/ exp(lGAMA)) + exp(lKS0)*log(DOSE)*exp(-exp(lGAMA)*T))),
               nl = TRUE,
               lKG ~ 1 + (1 | ID),
               lKS0 ~ 1 + (1 | ID),
               lGAMA ~ 1 + (1 | ID)
)
analytic

priors<- prior(normal(0,1),nlpar="lKG")+
  prior(normal(0,1),nlpar="lKS0")+
  prior(normal(0,1),nlpar="lGAMA")+
  prior(normal(0,1),class="sd",nlpar="lKG")+
  prior(normal(0,1),class="sd",nlpar="lKS0")+
  prior(normal(0,1),class="sd",nlpar="lGAMA")
priors

fit1<-brm(analytic,data=new_data,prior=priors,family=gaussian(),warmup = 1000, iter = 3000, 
          cores = 4, chains = 4, 
          seed = 123)
summary(fit1)
plot(fit1)
plot(conditional_effects(fit1),points=TRUE)

pp_check(fit1)















#-------------------ODE form: Simulation---------------------------------------
#------------------------------------------------------------------------------

# First try 


# MODEL: parameters, initial values, ODE
set.seed(123)

ODEcode<-'

$PARAM  TVKG = 0.001, TVKS0=0.2, TVGAMA=0.9, DOSE=10

$INIT
RESP=70

$MAIN


double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));



$OMEGA 0.005 0.03 0.003


$ODE

double KS=KS0 *log(DOSE)*exp( -GAMA * SOLVERTIME);

dxdt_RESP = KG * RESP -   KS* RESP ;

$CAPTURE KS GAMA 
'

set.seed(123)
ODEmod<-mcode("a1",ODEcode)
ODEmod




#' Time matrix
#'----------------------------------------------------
set.seed(123)
func_weeks_ODE<-function(nTAm, skel){
  rTAdys<-round(skel,2)
  #' Jitter TA days by 'amount' days
  TAweek<-round(rTAdys, 3)              #' Round to 3 digit years
  return(data.frame(TAnum=0:nTA, TAweek))
} 

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5
nTA<-9
skel<-seq(0,1,1/9)

set.seed(123)
indf1<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~func_weeks_ODE(nTAm = nTA,skel = skel), .id="ID")




#' Design matrix: Flat dose
#'----------------------------------------------------

#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
#' and to each ID corresponds a weight

set.seed(123)
doses<-c(10, 10, 10, 10, 10)

flt_ODE<-indf1 %>%
  mutate(TIME=TAweek, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=nTA+1))



set.seed(123)
snd.flt_ODE<-mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="TIME,KS,DOSE") %>% as.data.frame()

mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="TIME") %>%  plot(RESP+KS~TIME)




