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
  
  real dose[N];// covariates dose 
  
  int<lower=1> nTime[J];// number of observations per subject
  int<lower=1> T;// number of timepoints
  real ts[T];//timepoints for which ODE solutions are requested
  
  
}

transformed data{
  real x_r[1] ;
  int x_i[0];
}


// The parameters accepted by the model. 
// we want to parametrize the  inter-individual variability
parameters {
  
  // population parameter
  real<lower=0> tvkg;
  real<lower=0> tvks0;
  real<lower=0> tvgamma;
  real<lower=0> tvalpha;
  real<lower=0> tvbase;
  real<lower=0> sigma;
  
  real<lower=0> eta_kg[J];
  real<lower=0> eta_ks0[J];
  real<lower=0> eta_gamma[J];
  real<lower=0> eta_alpha[J];
  real<lower=0> eta_base[J];
  
  
  //omega
  real<lower=0> omega_kg;
  real<lower=0> omega_ks0;
  real<lower=0> omega_gamma;
  real<lower=0> omega_alpha;
  real<lower=0> omega_base;
  

  // thetas
  
  
}


transformed parameters{
  
  real<lower=0> kg[J];
  real<lower=0> ks0[J];
  real<lower=0> gamma[J];
  real<lower=0> alpha[J];
  real<lower=0> base[J];
  
  
  for(j in 1:J){
  kg[j]=exp(log(tvkg)+eta_kg[j]);
  ks0[j] =exp(log(tvks0)+eta_ks0[j]);
  gamma[j] =exp(log(tvgamma)+eta_gamma[j]);
  alpha[j] =exp(log(tvalpha)+eta_alpha[j]);
  base[j] =exp(log(tvbase)+eta_base[j]);
}
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  real y_hat[T,1];
  int index;
 
  //prior distributiona
  
  
  tvkg~lognormal(3.5,1);
  tvks0~beta(1,100);
  tvgamma~beta(1,100);
  tvalpha~beta(1,100);
  tvbase~beta(2,4);
  
  sigma~lognormal(0,1);
  
  
  
  
  omega_kg~lognormal(0,1);
  omega_ks0~lognormal(0,1);
  omega_gamma~lognormal(0,1);
  omega_alpha~lognormal(0,1);
  omega_base~lognormal(0,1);
  
  
  
  index=0;
  
  for(j in 1:J){
    real y0[1];
    real theta[4];
    
    eta_kg[j] ~ lognormal(0,omega_kg);
    eta_ks0[j] ~ lognormal(0,omega_ks0);
    eta_gamma[j] ~ lognormal(0,omega_gamma);
    eta_alpha[j] ~ lognormal(0,omega_alpha);
    eta_base[j] ~ lognormal(0,omega_base);
    
    //initial conditions
    y0[1]=base[j];
    theta[1]=kg[j];
    theta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    
    
    
  
    //likelihood
    for(t in 1:nTime[1]){
      
    index=index+1;
    
    y_hat = integrate_ode_rk45(sld, y0, t0,  ts,  theta,  x_r,  x_i);
    print(t0);
    
   // y_pred[index]= y_hat;
    y[index] ~ normal(y_hat[index],sigma);
      
    }
  }
}







