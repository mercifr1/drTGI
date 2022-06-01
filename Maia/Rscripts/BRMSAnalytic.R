source("C:/Users/muresai1/Desktop/drTGI/Maia/Rscripts/SimPRED.R")


#'Fit with brms 
#'-----------------------


#'Prepare the data for fit
#'
#'
#'Flat Dose:
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
new_data<-snd.flt%>%
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMMA=log(GAMMA),lBAS=log(BASE))%>%subset(select=-c(KS0,KG,GAMMA,BASE))
names(new_data)

#' Define the model structure and parameters < random effects>
analytic <- bf(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMA) ) *( 1 - exp(exp(lGAMA)*TIME)))  ,
               nl = TRUE, 
               lBAS ~ 1+ (1|ID),
               lKG ~ 1 + (1 | ID),
               lKS0 ~ 1 + (1 | ID),
               lGAMA ~ 1 + (1 | ID))
analytic


#'** Define the priors**
#'
#'
#'
# check which  priors
getpr1<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = new_data, family = gaussian())


getpr2<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = new_data, family = lognormal())

getpr3<-get_prior(RESP ~ exp(lBAS) *  exp( exp(lKG) * TIME - ( exp(lKS0) *log(DOSE)/ exp(lGAMMA) ) *( 1 - exp(exp(lGAMMA)*TIME)))  ,
                  data = new_data, family = student())

# # first try: it is not converging:max(Rhat) is  > 3
# 
# prior1<-set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lBAS")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lKG")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lKS0")+
#   set_prior("student_t(3, 67.2, 10)",class="b",nlpar="lGAMA")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lBAS")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lKG")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lKS0")+
#   set_prior("student_t(3, 0, 10)",class="sd",nlpar="lGAMA")
# prior1
# 
# m01<-brm(analytic, data=new_data, prior=prior1, family=gaussian(),
#          
#          iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
# plot(m01)



# second try:It converge: Rhat=1 

#' Define the model structure and parameters < random effects>
analytic2 <- bf(RESP ~ exp(lBAS) *  exp( (exp(lKG)/100) * TIME - ( (exp(lKS0)/100) *log(DOSE)/ (exp(lGAMMA)/100) ) *( 1 - exp((exp(lGAMMA)/100)*TIME)))  ,
                nl = TRUE, 
                lBAS ~ 1+ (1|ID),
                lKG ~ 1 + (1 | ID),
                lKS0 ~ 1 + (1 | ID),
                lGAMMA ~ 1 + (1 | ID))
analytic2

prior2<-set_prior("normal(4,1)",class="b",nlpar="lBAS")+
  set_prior("normal(0,1)",class="b",nlpar="lKG")+
  set_prior("normal(0,1)",class="b",nlpar="lKS0")+
  set_prior("normal(0,1)",class="b",nlpar="lGAMMA")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lBAS")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKG")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lKS0")+
  set_prior("student_t(3, 0, 2.5)",class="sd",nlpar="lGAMMA")
prior2

m02<-brm(analytic2, data=new_data, prior=prior2, family=lognormal(link = "identity"),
         
         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m02)
plot(m02)

pairs(m02)

#' 2 Dose levels:
#'-----------------------------------------------------------------------------
#' log-transformation of the parameters
#' 


new_data_2<-snd.tv%>%
  mutate(lKG=log(KG),lKS0=log(KS0),lGAMA=log(GAMMA),lBAS=log(BASE))%>%subset(select=-c(KS0,KG,GAMMA,BASE))
names(new_data_2)
head(new_data_2)

m03<-brm(analytic2, data=new_data_2, prior=prior2, family=lognormal(link = "identity"),
         
         iter=2000, chains=4, warmup=1000, seed=1234, cores=4,control = list(adapt_delta = .9))
summary(m03)
plot(m03)




