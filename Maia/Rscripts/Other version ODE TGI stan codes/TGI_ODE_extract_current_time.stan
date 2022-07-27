

functions{
  real[] sld(real t, real[] y, real[] theta, real[] x_r, int[] x_i){
    
    real dydt[1];
    //real dose=x_r[1];
  
    // define the ODE
     
    dydt[1]=theta[1]*y[1]-theta[2]*exp(-theta[3]*t)*theta[4]*y[1];
    return dydt;
  }
  
  real[] extract_ts(real[] ts, int j, int index, int[] nTime) {
    real return_obj[nTime[j]];
    
    return_obj = ts[index:(index + nTime[j] -1)];
    
    return return_obj;
  }
}




// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;//total number of observations
  int<lower=1> J;// number of individuals
  int<lower=1> nTime[J];// number of observations per subject 
  int<lower=1> uniqueT;// number of timepoints all unique time points
  real ts[uniqueT];//timepoints for which ODE solutions are requested
  int<lower=0> timeInd[N]; // numeric time indicator
  real y[N];//response vector
  real t0;// starting time
  
  //real dose[N];// covariates dose 
}

transformed data{
  real x_r[0] ;
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
  
  real eta_kg[J];
  real eta_ks0[J];
  real eta_gamma[J];
  real eta_alpha[J];
  real eta_base[J];
  
  
  //omega
  real<lower=0> omega_kg;
  real<lower=0> omega_ks0;
  real<lower=0> omega_gamma;
  real<lower=0> omega_alpha;
  real<lower=0> omega_base;
  
  
}

transformed parameters{
  real<lower=0> kg[J];
  real<lower=0> ks0[J];
  real<lower=0> gamma[J];
  real<lower=0> alpha[J];
  real<lower=0> base[J];
  
  for( j in 1:J){
   kg[j] = exp(tvkg+eta_kg[j]*omega_kg);
   ks0[j] =exp(tvks0+eta_ks0[j]*omega_ks0);
   gamma[j] =exp(tvgamma+eta_gamma[j]*omega_gamma);
   alpha[j]= exp(tvalpha+eta_alpha[j]*omega_alpha);
   base[j] = exp(log(tvbase)+eta_base[j]*omega_base);
  }
}
 


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  real y_hat[uniqueT,1];
  int index;
 
  //prior distributiona
  
  tvkg~normal(0.5,0.5);
  tvks0~normal(0.5,0.5);
  tvgamma~normal(0.5,0.5);
  
  tvalpha~normal(0,2);
  
  tvbase~normal(80,10);

  sigma~normal(0,1);

  omega_kg~normal(0,1);
  omega_ks0~normal(0,1);
  omega_gamma~normal(0,1);
  omega_alpha~normal(0,1);
  omega_base~normal(0.5,0.5);



  eta_kg ~ normal(0,1);
  eta_ks0 ~ normal(0,1);
  eta_gamma ~ normal(0,1);
  eta_alpha ~ normal(0,1);
  eta_base ~ normal(0,1);

  
  index=0;
  
  for(j in 1:J){
    real y0[1];
    real theta[4];
    real curr_ts[nTime[j]];
    
    
    //initial conditions
    y0[1]=base[j];
    theta[1]=kg[j];
    theta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    
    
   
    
    curr_ts = extract_ts(ts, j, index + 1, nTime);
    y_hat[(index + 1 ):(index + nTime[j] )] = integrate_ode_rk45(sld, y0, t0, curr_ts, theta, x_r, x_i);
   
  
    //likelihood
    for(t in 1:nTime[j]){
      
    index = index+1;
    
   // y_pred[index]= y_hat;
    y[index] ~ normal(y_hat[timeInd[index],1],sigma);
      
    }
  }
}





