source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimAnalyticalwithTIME.R")
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimODEwithSOLVERTIME.R")


#' Design matrix: 2 dose levels at a certain time point for all patients
#' ---------------------------------------------------

stepup0<-flt%>%mutate(DOSE=ifelse(flt$tnum==5,5,10))
head(stepup0,n=10)


#' **2. For the 2 level dose**
#' 
#' 
#' 
#' For PRED with SIGMA:
set.seed(123)
snd.stepup0<-mrgsim_d(model_pred_sigma, as.data.frame(stepup0), carry.out="DOSE") %>% 
  as.data.frame()
  head(snd.stepup0,n=20)


#' 
#' 
#' For ODE with SIGMA: 
#' 
set.seed(123)
snd.stepup0_ODE<-mrgsim_d(model_ode_sigma, as.data.frame(stepup0),carry.out="DOSE") %>% 
  as.data.frame()

head(snd.stepup0_ODE,n=20)




#' Plot both simulations

set.seed(123)
mrgsim_d(model_pred_sigma, as.data.frame(stepup0), carry.out="DOSE")  %>% plot(RESP~time)
set.seed(123)
mrgsim_d(model_ode_sigma, as.data.frame(stepup0), carry.out="DOSE")  %>% plot(DV~time)
#------------------------------------------------------------------------------------------------

