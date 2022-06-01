source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODE.R")

#'Scenario 1 : balanced design matrix

#' -   dose-cohorts with same number of individual/cohort:
#'       - 5 cohorts
#'       - 5 individual/cohort
#'      - 25 individual in total

#' -  each cohort have a different flat/constant dose=> 5 dose levels 
#' -  all individuals:  same number of observations 



#'Define the structure of the model
#' 

##'  with ODE

ode<-'
$PARAM  TVKG =0.6, TVKS0=0.4, TVGAMMA=0.8, TVBASE =70, DOSE=10

$OMEGA  0.005 0.03 0.003 0.01
$SIGMA  0.05

$CMT RESP

$MAIN
double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMMA = TVGAMMA * exp(ETA(3));
double BASE = TVBASE*exp(ETA(4));
RESP_0 = BASE;


$ODE
double KS=KS0 *exp( -GAMMA * SOLVERTIME);
dxdt_RESP = KG *RESP -   KS*log(DOSE)*RESP;


$CAPTURE 
  RESP_0 TVKG TVKS0 EPS(1) \\ add the TVKG and TVKS0 to 
'
#' Compile the model

set.seed(123)
model_ode<-mcode("ODE",ode)

#' Construction of the Time matrix
#'----------------------------------------------------
nInd<-25 #number of individual

# build a function to define the time 
funtime<-function(nTAm,skel){
  tweek<-jitter(skel,amount = 7)
  ty<-round(tweek/365,3)
  return(data.frame(TAnum=0:nTA,ty))
}

nTA<-9 # number of observations
skel_6w<-7*c(0,6,12,18,24,30,36,42,48,52) #' dose administration at every 6 weeks
skel_1<-7*c(0,6,12,18,24,30,36,42,48,52,58,64,70)
                  
 
indf<-tibble(ID=1:nInd) %>%
  split(.$ID) %>%
  map_dfr(., ~funtime(nTA=nTA,skel=skel_6w), .id="ID")
head(indf,n=15)


#' Design matrix: 5 individual per each cohort with Constant dose over time
#'----------------------------------------------------


#' **To each ID corresponds a different dose e.g. 5 ID = 5 doses**

doses<-10*c(1/6, 1/2, 1, 2, 3) 


# Design matrix Scenario 1
flt.S1<-indf %>%
  mutate(TIME=ty,
         ID=as.numeric(as.character(ID)),
         DOSE=rep(round(doses,2), each=(max(indf["TAnum"])+1)*5),
         TVKS0=ifelse(DOSE==1.67,0.8,ifelse(DOSE==30,0.3,0.4)),
         TVKG=ifelse(DOSE==1.67,0.3,ifelse(DOSE==30,0.8,0.6)),
         EVID=2,
         CMT=1)%>%
  subset(select=-c(ty))

head(flt.S1,n=15)


#'Simulation
set.seed(123)
ode.flt2<-mrgsim_d(model_ode,as.data.frame(flt.S1),carry.out="DOSE",end=-1,nocb=F)
  
ode.flt2<-ode.flt2%>%mutate(DV=RESP+EPS_1) #additive error 
head(ode.flt2,n=20)

#' plot
p1<-ggplot(ode.flt2, aes(TIME,DV,group=ID,color=factor(DOSE)))+
  geom_line(size=0.8)+
  theme_bw()+
  xlab("time")+
  ylab("SLD")+
  guides(fill=guide_legend(title="DOSE"))
p1




