###############################################################################



#     Check the non-centered parametrization for the parameters for Stan



################################################################################
# Hierarchical models can be made much more efficient by shifting the data’s correlation 
# with the parameters to the hyperparameters.
# see https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html

#' *  baseline y0*

# prior
tvbase<-rnorm(100,80,10)
# sd of random effect
omega_base<-abs(rnorm(100,0,1))

#random effects
eta_base<-rnorm(100,0,1)

# non-cenntered parametrization 
y0<-exp(tvbase+eta_base*omega_base)
y0
hist(y0)
# with this expression we force y0 to pe positive but one can also conclude from
# the histograme 
# from fitting the model in stan with this parametrization we conclude that MCMC diverge


# better option for y0 expression
y0<-exp(log(tvbase)+eta_base*omega_base)
y0
hist(y0)
# this log transformation is working much better
# one reason can be, because the prior tvbase follows a normal distribution with mean 80 and sd 10.


#' *  growth rate kg*
tvkg<-rnorm(100,0.5,0.5)
omega_kg<-abs(rnorm(100,0,1))
eta_kg<-rnorm(rnorm(100,0,1))
kg<-exp(tvkg+eta_kg*omega_kg)
hist(kg)
# in this case, the parametrization and the chosen distribution seem feasible 
# and giving results around the desired value for kg


#' *  shrinkage rate ks0*
tvks0<-rnorm(100,0.5,0.5)
omega_ks0<-abs(rnorm(100,0,1))
eta_ks0<-rnorm(rnorm(100,0,1))
ks0<-exp(tvks0+eta_ks0*omega_ks0)
hist(ks0)

#' *  resistance parameter gamma*
tvgamma<-rnorm(100,0.5,0.5)
omega_gamma<-abs(rnorm(100,0,1))
eta_gamma<-rnorm(rnorm(100,0,1))
gamma<-exp(tvgamma+eta_gamma*omega_gamma)
hist(gamma)


#' * strength parameter between dose and ks: alpha  *
tvalpha<-rnorm(100,0,2)
omega_alpha<-abs(rnorm(100,0,1))
eta_alpha<-rnorm(rnorm(100,0,1))
alpha<-exp(tvalpha+eta_alpha*omega_alpha)
hist(alpha)