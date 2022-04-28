#Libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(purrr)



library(brms)


#'-------------------Analytically form: Simulation------------------------------
#'------------------------------------------------------------------------------

#' First try
#' --------------------------------------------------------------------


set.seed(123)

#' Implement the structure of the model

code <- "

$PARAM    TVKG =0.06, TVKS0=0.4, TVGAMA=0.8, TVBASE =70, DOSE=10, T=0


$OMEGA  0.005 0.03 0.003 0.01

$SIGMA 0.05

$PRED

double KG = TVKG * exp(ETA(1));
double KS0 = TVKS0 * exp(ETA(2));
double GAMA = TVGAMA * exp(ETA(3));

double BASE=TVBASE * exp(ETA(4));
capture RESP_0=BASE;

double KS =  KS0 * log(DOSE) * exp(-GAMA * T);

capture RESP = RESP_0*(exp((KG*T)-( KS0 * log(DOSE)/GAMA)+KS))+EPS(1);


$CAPTURE  RESP_0 BASE KS0 KS KG GAMA TIME 
"


#' Compile the model code
#' ------------------------------------------------------

set.seed(123)
mod <- mcode("maia", code)
mod


#' Time matrix
#'----------------------------------------------------
#'
#'
set.seed(123)

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5
funtime<-function(n,end,delta){
  ty<-seq(0,end,delta)
  return(data.frame(tnum=0:n,ty))
}

indf2<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~funtime(n=20,end=1,delta=0.05) , .id="ID")
indf2


#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
doses<-c(10, 10, 10, 10, 10)

flt2<-indf2 %>%
  mutate(T=ty, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=20+1))%>%subset(select=-c(ty))
flt2



#' Simulation 
#' ---------------------------------------------------
set.seed(123)
snd.flt2<-mrgsim_d(mod, as.data.frame(flt2), carry.out="KS,KS0,DOSE,RESP_0,BASE,T") %>% 
  as.data.frame()%>%subset(select=-c(time,TIME))
head(snd.flt2,n=10)



ggplot(snd.flt2, aes(T,RESP,group=ID)) + geom_line()




#########################################
set.seed(123)

#' **Nb of Ind must corresponds to nb of doses in this program**
nInd<-5
funtime<-function(end,delta){
  ty<-seq(0,end,delta)
  p<-length(ty)
  return(data.frame(num=0:(p-1),ty))
}
funtime(24,0.1)


indf3<-tibble(ID=1:nInd) %>% 
  split(.$ID) %>% 
  map_dfr(., ~funtime(end=24,delta=1) , .id="ID")
indf3


#' To each ID corresponds a different dose e.g. 5 ID = 5 doses
doses<-c(10, 10, 10, 10, 10)

flt3<-indf3 %>%
  mutate(T=ty, ID=as.numeric(as.character(ID)),
         DOSE=rep(doses, each=24+1))%>%subset(select=-c(ty))
flt3



#' Simulation 
#' ---------------------------------------------------
set.seed(123)
snd.flt3<-mrgsim_d(mod, as.data.frame(flt3), carry.out="KS,KS0,DOSE,RESP_0,BASE,T") %>% 
  as.data.frame()%>%subset(select=-c(time,TIME))
head(snd.flt3,n=10)

mrgsim_d(mod, as.data.frame(flt3), carry.out="KS,KS0,DOSE,RESP_0,BASE,T") %>% plot(RESP~T)

ggplot(snd.flt3, aes(T,RESP,group=ID)) + geom_line()



