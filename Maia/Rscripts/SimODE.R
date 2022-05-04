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
KS  RESP_0 BASE KG GAMA
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
mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="DOSE,T",zero_re(sigma)) %>%  plot(RESP~T)

out<-mrgsim_d(ODEmod, as.data.frame(flt_ODE), carry.out="DOSE,T") %>%as.data.frame()
out$ID<-as.factor(out$ID)

ggplot(out,aes(T,RESP,color=ID))+
  geom_line()+
  scale_y_continuous(breaks = seq(50,100,5))

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

