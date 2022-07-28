################################################################################

# Fit TGI model:  a mixed-effect model defined with one ODE 

# Solving an  hierarchical ODE for time arrays of different size with NLMIXR

# nlmixr R package for fitting general dynamic models,PK models and PKPD models in particular

# see :https://nlmixrdevelopment.github.io/nlmixr/


# 31-07-2022

# Maia Muresan

################################################################################

source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/Simulation_SLD_complete.R")


# libraries
#install.packages("nlmixr")
library(nlmixr)
#install.package("broom.mixed")
library(broom.mixed)


# prepare data set for nlmixr

nlmixr_data<- design_flat%>% # design_flat data set from Simulation_SLD_complete
  filter(REP==1)%>%
  mutate(id=ID, time=ifelse(TIME< 0.1, 0 , TIME), dv=DV, EVID = 0, FLAG = 2, cmt = 1, dvid = "cp")%>%
  select( id, time, dv, DOSE, EVID, cmt, dvid)

nlmixr_data$FLAG<- as.integer(nlmixr_data$FLAG)
nlmixr_data$EVID<- as.integer(nlmixr_data$EVID)
nlmixr_data$id<- as.integer(nlmixr_data$id)



#define the model 

tgiModel<-function(){ #Claret-like model for describing SLD
  ini({
    tKG <- log(0.6) # logKG
    tKS0 <- log(0.2) # logKS0
    tGAMMA <- log(0.8)# logGAMMA
    tALPHA <- log(2) #logALPHA
    tBASE<- log(70)# logIBASE
    
    eta.KG ~ 0.005
    eta.KS0 ~ 0.03
    eta.GAMMA ~ 0.003
    eta.ALPHA ~ 0.4
    eta.BASE ~ 0.01
    
    add.err<-20 #SD of additive variability, Sd of additive error
  })
  model({
    KG=exp(tKG + eta.KG) # tumor growth rate constant
    KS0=exp(tKS0 + eta.KS0) # shrinkage  rate constant
    GAMMA=exp(tGAMMA + eta.GAMMA) 
    ALPHA=exp(tALPHA + eta.ALPHA)
    BASE=exp(tBASE + eta.BASE) # baseline sld (mm)
    
    SLD(0)=BASE*1000 #baseline sld values
    d/dt(SLD)=KG * SLD - (KS0/1000) * exp((GAMMA/1000)*time) * (ALPHA/1000) * log(1+DOSE) * SLD # Model for dSLD_dt
    cp = SLD
    cp ~ add(add.err)
  })
  
}


fit_nlmixr<-nlmixr(tgiModel, data=nlmixr_data, est="saem", control=list(print=0),
                   table=list(cwres=TRUE,npde=TRUE))
print(fit_nlmixr)
fit_nlmixr$covMethod
fit_nlmixr<- fit_nlmixr %>% addCwres()%>% addNpde()

plot(fit_nlmixr)


glance(fit_nlmixr)
setOfv(fit_nlmixr,"gauss3_1.6")
tidy(fit_nlmixr)

## transformation applied on every parameter
tidy(fit_nlmixr,exponentiate=TRUE)
## transformation applied to log-scaled population parameters


tidy(fit_nlmixr, effects="fixed")
tidy(fit_nlmixr, effects="ran_pars")

head(tidy(fit_nlmixr, effects="ran_vals"))
head(tidy(fit_nlmixr, effects="ran_coef"))