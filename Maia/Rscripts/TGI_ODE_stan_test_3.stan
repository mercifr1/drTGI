//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled 
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions{
  real[] sld(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    
    real dydt[1];
    real dose=x_r[1];
  
    // define the ODE
     
    dydt[1]=theta[1]*y[1]-theta[2]*exp(-theta[3]*t)*theta[4]*log(1+dose)*y[1];
    return dydt;
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;//total number of observations
  int<lower=1> J;// number of individuals
  real y[N];//response vector
  real t0;// starting time
  
  real dose;// covariates dose 
  int<lower=1> timeInd[N]; // numeric time indicator
  int<lower=1> nTime[J];// number of observations per subject
  int<lower=1> T;// number of timepoints
  real ts[T];//timepoints for which ODE solutions are requested
  
  
}

transformed data{
  real x_r[1] = {dose}  ;
  int x_i[0];
}


// The parameters accepted by the model. 
// we want to parametrize the  inter-individual variability
parameters {
  
  // population
  real<lower=0> tvkg;
  real<lower=0> tvks0;
  real<lower=0> tvgamma;
  real<lower=0> tvalpha;
  real<lower=0> tvbase;
  real<lower=0> sigma;
  
  //etas
  // real<lower=0> eta_kg;
  // real<lower=0> eta_ks0;
  // real<lower=0> eta_gamma;
  // real<lower=0> eta_alpha;
  // real<lower=0> eta_base;
  
  //omega
  real<lower=0> omega_kg;
  real<lower=0> omega_ks0;
  real<lower=0> omega_gamma;
  real<lower=0> omega_alpha;
  real<lower=0> omega_base;
  
  
  
  
  // thetas
  vector[J] kg;
  vector[J] ks0;
  vector[J] gamma;
  vector[J] alpha;
  vector[J] base;
  
}


transformed parameters{
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real y_pred[N];
  int index;
  
  //prior distributiona
  tvkg~lognormal(0.2,0.05);
  tvks0~lognormal(0.2,0.05);
  tvgamma~lognormal(0.2,0.05);
  tvalpha~lognormal(0.2,0.05);
  tvbase~lognormal(0.2,0.05);
  
  omega_kg~inv_gamma(2,0.1);
  omega_ks0~inv_gamma(2,0.1);
  omega_gamma~inv_gamma(2,0.1);
  omega_alpha~inv_gamma(2,0.1);
  omega_base~inv_gamma(2,0.1);
  
  
  index=0;
  for(j in 1:J){
    real y0[1];
    real theta[4];
    
    
    //initial conditions
    y0[1]=base[j];
    
    theta[1]=kg[j];
    theta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    
    
    
    //likelihood
    for(t in 1:nTime[J]){
    index=index+1;
    
  real y_hat[1,1] =  integrate_ode_rk45(sld, y0, t0,  ts,  theta,  x_r,  x_i);
    
     y_pred[index]= y_hat[1,1];
      
      
    y[index] ~ normal(y_pred[timeInd[index]],sigma);
    }
    kg[j] ~  normal(tvkg,omega_kg);
    ks0[j] ~  normal(tvks0,omega_ks0);
    gamma[j] ~  normal(tvgamma,omega_gamma);
    alpha[j] ~  normal(tvalpha,omega_alpha);
    base ~  normal(tvbase ,omega_base);
    
    
    
    
  }
     
     
}

