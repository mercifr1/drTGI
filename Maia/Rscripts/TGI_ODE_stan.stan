

functions {
  real sld(real t,
  real y,
  real kg,
  real ks0,
  real gamma,
  real base,
  real alpha,
  real x_r,
  int x_i
  ){
    
    
    real dose = x_r;
    real dy_dt=kg*y-ks0*exp(-t*gamma)*alpha*log(1+dose)*y;
    return dy_dt;
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
  real<lower=0> kg;
  real<lower=0> ks0;
  real<lower=0> gamma;
  real<lower=0> alpha;
  real<lower=0> base;
  real<lower=0> sigma;
}

transformed parameters{
  
}

model {
  real y_hat[T,1];
  int index;
  //prior distributions
  kg~lognormal();
  ks0~lognormal();
  gamma~lognormal();
  base~lognormal();
  alpha~lognormal();
  
  index=0;
  for(j in 1:J){
    real y0;
  }
}


