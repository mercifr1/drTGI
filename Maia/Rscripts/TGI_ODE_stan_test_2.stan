

functions {
  real sld(real t,
  real y,
  real kg,
  real ks0,
  real gamma,
  real base,
  real alpha,
  real dose,
  real x_r,
  int x_i
  ){
    
    real dy_dt=kg*y-ks0*exp(-t*gamma)*alpha*log(1+dose)*y;
    return dy_dt;
  }
}

data {
  int<lower=1> N;//total number of observations
  int<lower=1> J;//total nomber of individual
  real y[N]; //response vector
  real t0; // starting time
  int<lower=1> T; // number of timepoints (all unique timepoints)
  int<lower=1> nTime[J]; // number of observations per individual
  int<lower=1> timeID[N]; //numeric time indicator
  real ts[T];//timepoints for which ODE solutions are requested
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  vector[J] kg;
  vector[J] ks0;
  vector[J] gamma;
  vector[J] alpha;
  vector[J] base;
  vector[N] sigma;
  
}

transformed parameters{
  
}

model {
  real y_hat[T,1];
  int index;
  //prior distributions
  kg[J]~lognormal(0.6,0.005);
  ks0[J]~lognormal(0.2,0.03);
  gamma[J]~lognormal(0.8,0.003);
  base[J]~lognormal(70,0.01);
  alpha[J]~lognormal(2,0.04);
  
  index=0;
  for(j in 1:J){
    real y0;
    real theta[4];
    //initial conditions
    y0=base[j];
    theta[1]=kg[j];
    theta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    
    y_hat=integrate_ode_rk45(sld, y0, t0, ts, theta, x_r, x_i);
    
    //likelohood
    
    
  }
}


