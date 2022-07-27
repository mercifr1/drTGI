###############################################################################


# TGI: Fit the Analytically form of the ODE with brms package

# BRMS: Bayesian regression models using Stan see:https://paul-buerkner.github.io/brms/

# 31-07-2022

# Maia Muresan

################################################################################



# source the analytical form simulation from SimPRED.R
source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimPRED.R")

# not used here: in case we want to take the complete simulations 
#source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Simulation_SLD_complete.R")

#'Fit with brms 
#'-----------------------
library(brms)

#'Prepare the data for fit
#'
#'
#'Flat Dose:
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
#' 

brms_data<-snd.flt%>%# snd.flt data set is from SimPRED.R 
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMMA=log(GAMMA),lBAS=log(BASE))%>%
  subset(select=-c(KS0,KG,GAMMA,BASE))
names(brms_data)

#' Define the model structure and parameters < random effects>
#' # It converge: Rhat=1 


analytic <- bf(RESP ~ exp(lBAS) *  exp( (exp(lKG)/100) * TIME - ( (exp(lKS0)/100) *log(DOSE)/ (exp(lGAMMA)/100) ) *( 1 - exp((exp(lGAMMA)/100)*TIME)))  ,
                nl = TRUE, 
                lBAS ~ 1+ (1|ID),
                lKG ~ 1 + (1 | ID),
                lKS0 ~ 1 + (1 | ID),
                lGAMMA ~ 1 + (1 | ID))
analytic

prior<-set_prior("normal(4,1)",class="b",nlpar="lBAS")+
  set_prior("normal(0,1)",class="b",nlpar="lKG")+
  set_prior("normal(0,1)",class="b",nlpar="lKS0")+
  set_prior("normal(0,1)",class="b",nlpar="lGAMMA")+
  
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lBAS")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKG")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKS0")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lGAMMA")


prior

m01<-brm(analytic, data=brms_data, prior=prior, family=lognormal(link = "identity"),
         
         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m01)

stancode(m01)
plot_brms_PRED_flt<-plot(m01)

pairs(m01)

#'** Define the priors**
#'
#'
#'
# check which  priors
getpr1<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = brms_data, family = gaussian())


getpr2<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = brms_data, family = lognormal())

getpr3<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = brms_data, family = student())







#' time varying dose: 2 level dose:  10 and 5 
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
#' 


brms_data_2<-snd.tv%>% # snv.tv data set is from SimPred.R
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMA=log(GAMMA),lBAS=log(BASE))%>%subset(select=-c(KS0,KG,GAMMA,BASE))
names(new_data_2)
head(new_data_2)

m02<-brm(analytic, data=brms_data_2, prior=prior, family=lognormal(link = "identity"),
         
         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m02)
plot_brms_PRED_tv<-plot(m02)


make_stancode(analytic,data=brms_data_2,family = gaussian())
