//                       Fitting TGI model in Stan 

// Remark: Cmdstan 2.24 introduce a new ODE interface see: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
// The Stan code above is implemented in the old solver inference 




// First step: Coding the ODE System Function
// see https://mc-stan.org/docs/2_25/stan-users-guide/coding-the-ode-system-function.html


functions{
  // Important: The types in the ODE system function are strict and require a specific signature
  
  
  // t: first argument is the  single time point passed as <real> at which we want to evaluate the derivative.
  // y: second argument is the state i.e the scalar function value that in our case will be the SLD values, passed as an  real array or vector , 
  //    and the return value from the system function are the current time derivatives of the state defined as an array
  // theta: parameter vector (4 parameters).
  // x_r : real covariates from the data : in this case the < DOSE >
  // x_i:  integer data : not used.
  
  // The function returns a vector of time derivatives of the system state at time t, state y and parameter theta
  
  // Function signature needs to be conformant with what integrate_ode_* accepts.
  real[] dSLD_dt(real t, // time
                 real[] y, // state
                 real[] theta, // parameters ( 4 parameters)
                 real[] dose, // real data : covariate dose in our case
                 int[] x_i){  // integer data
    
    real dydt[1];
    
  
    // define the ODE 
    // dydt=kg*y-ks0*exp(-gamma *t)*alpha*log(1+dose)*y
    dydt[1]=theta[1]*y[1]-theta[2]*exp(-theta[3]*t)*theta[4]*log(1+dose[1])*y[1];
    return dydt;
  }
}




// The input data 
data {
  int<lower=0> N;   //total number of observations 
  int<lower=0> J;   // total  number of individuals
  int<lower=0> nTime[J]; // number of observations for each  subject 
  int<lower=0> uniqueT; // total number of timepoints: all unique time points for all subjects
  real ts[uniqueT];   //timepoints for which ODE solutions are requested
  real y[N];         //response vector
  real t0;          // starting time
  
  real dose[uniqueT];// covariates dose 
}

transformed data{
  int x_i[0]; // integer data is of length 0 because is not used here
  
  //compute index showing from which time position the solution of ODE system should be derived for each subject.
  // useful later in the code for the integrator_ode*specific : 
  // - vector of time poisnt <ts[..]> from data block contains different time observations for each subject,
  //   and the integrator must conatin a strictly increasing vector 
 
  int t_index[J]; 
  t_index[1]=1;
  for (j in 2:J){
    t_index[j]=t_index[j-1]+nTime[j-1];
  }
}


// The parameters accepted by the model. 
// we want to parametrize the  inter-individual variability
parameters {
  
  // population parameter fixed-effects (population mean of parameters) 
  real<lower=0> tvkg;
  real<lower=0> tvks0;
  real<lower=0> tvgamma;
  real<lower=0> tvalpha;
  real<lower=10> tvbase;
  
  // sd of the residual error, lower=0 for half-normal
  real<lower=0> sigma;
  
  //omega: sd of the random-effects for  parameters), lower=0 for half-normal
  real<lower=0> omega_kg;
  real<lower=0> omega_ks0;
  real<lower=0> omega_gamma;
  real<lower=0> omega_alpha;
  real<lower=0> omega_base;
  
  
  //etas:  random effects of parametrers for each individual
  real eta_kg[J];
  real eta_ks0[J];
  real eta_gamma[J];
  real eta_alpha[J];
  real eta_base[J];
  
  
  
  
}

transformed parameters{
  
  real IPRED[uniqueT];// individual prediction

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

  
  
  for(j in 1:J){
    
    
    real y0[1];
    real theta[4];
    
    real sol[nTime[j],1];
    
    //initial conditions
    y0[1]=base[j];
    theta[1]=kg[j];
    theta[2]=ks0[j];
    theta[3]=gamma[j];
    theta[4]=alpha[j];
    
    // see https://mc-stan.org/docs/functions-reference/functions-old-ode-solver.html
    //  *vector segment( vector v, int i, int n)* : return the row vector consisting of the n elements of v starting from i
    // see https://mc-stan.org/docs/2_26/functions-reference/slicing-and-blocking-functions.html
    sol = integrate_ode_bdf(dSLD_dt, y0, t0, segment(ts,t_index[j] , nTime[j]), theta, dose, x_i);
    
  
    // individual predictions
    for(t in 1:nTime[j]){
    IPRED[t_index[j] + t-1]=sol[t,1];}
    
  
   }
}
 


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'IPRED'
// and standard deviation 'sigma'.


model {
  
  //prior distributions
  
 
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

  
  //Model for the data (likelihood)
  y ~ normal(IPRED,sigma);
      
}




