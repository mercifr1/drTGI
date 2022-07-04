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
functions{
  real[] sld(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    real dydt[1];
    dydt[1]=theta[1]*y[1]-theta[2]*exp(-theta[3]*t)*theta[4]*log(1+x_r[1])*y[1];
    return dydt;
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;//total number of observations
  int<lower=1> J;// number of individuals
  real y[N];//response vector
  real t0;// starting time
  //real dose[N];// dose for each observation
  int<lower=1> nTime[J];// number of observations per subject
  real ts;
  int<lower=1> T;// number of timepoints
  
}

transformed data{
  real x_r[1];
  
  real x_i[0];
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  // 
  real<lower=0> tvkg;
  real<lower=0> tvks0;
  real<lower=0> tvgamma;
  real<lower=0> tvalpha;
  real<lower=0> tvbase;
  real<lower=0> sigma;
  
  // etas
  real<lower=0> eta_kg;
  real<lower=0> eta_ks0;
  real<lower=0> eta_gamma;
  real<lower=0> eta_alpha;
  real<lower=0> eta_base;
  
  real kg[J];
  real ks0[J];
  real gamma[J];
  real alpha[J];
  real base[J];
  
}


transformed parameters{
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real y_hat[T];
  int index;
  
  //prior distributiona
  tvkg~lognormal(0.2,0.05);
  tvks0~lognormal(0.2,0.05);
  tvgamma~lognormal(0.2,0.05);
  tvalpha~lognormal(0.2,0.05);
  tvbase~lognormal(0.2,0.05);
  
  eta_kg~inv_gamma(2,0.1);
  eta_ks0~inv_gamma(2,0.1);
  eta_gamma~inv_gamma(2,0.1);
  eta_alpha~inv_gamma(2,0.1);
  eta_base~inv_gamma(2,0.1);
  
  
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
    
    y_hat=integrate_ode_rk45(sld,y0,t0,ts,theta,x_r,x_i);
    
    //likelihood
    for(t in 1:nTime[J]){
      index=index+1;
      y[index]~normal(y_hat[timeInd[index]],sigma);
    }
    kg[j]~normal(tvkg,eta_kg);
    ks0[j]~normal(tvks0,eta_ks0);
    gamma[j]~normal(tvgamma,eta_gamma);
    alpha[j]~normal(tvalpha,eta_alpha);
    base~normal(tvbase,eta_base);
    
    
    
    
  }
     
     
}

