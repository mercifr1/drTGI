//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real[] sld(real t,real[] y,real[] theta,real[] x_r,int[] x_i){
    real dydt;
    real dose = x_r;
    dydt=theta[1]*y-theta[2]*exp(-theta[3]*t)*theta[4]*log(1+dose)*y;
    return dydt;
  }
}

data{
  int<lower=1> N;//total number of observations
  int<lower=1> J;//total nomber of individual
  real y[N]; //response vector
  real t0; // starting time
  
  int<lower=1> T; // number of timepoints (all unique timepoints)
  int<lower=1> nTime[J]; // number of observations per individual
  int<lower=1> timeID[N]; //numeric time indicator
  real ts[T];//timepoints for which ODe solutions are requested
}

transformed data {
  real x_r[1];
  int x_i[0];
}

parameters {
  real<lower=0> kgm;
  real<lower=0> ks0m;
  real<lower=0> gammam;
  real<lower=0> alpham;
  real<lower=0> basem;
  real<lower=0> sigmam;
  real kg[J];
  real ks0[J];
  real gamma[J];
  real alpha[J];
  real base[J];
}

model {
  real y_hat[T,1];
  int index;
  //prior distributions
  kgm~lognormal();
  ks0m~lognormal();
  gammam~lognormal();
  basem~lognormal();
  alpham~lognormal();
  sigmam~lognormal();
  
  index=0;
  for(j in 1:J){
    real y0[1];
    real theta[5];
    //initial conditions
    y0=base[j];
    theta[1]=kg[j];
    tehta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    y_hat=integrate_dose_rk45(sld,y0,t0,ts,theta,x_r,x_i);
  }
}
