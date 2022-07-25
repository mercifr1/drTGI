

//

functions {
  //The types in the ODE system function are strict:
  // t: first argument is the  single time point passed as <real> at which we want to evaluate the derivative.
  // sld: second argument is the state i.e the scalar function value, passed as an <array>, and the return value
  //      from the system function are the current time derivatives of the state defined as an array
  // theta: parameter vector (4 parameters).
  // x_r : real covariates from the data : in this case the < DOSE >
  //x_i:  integer data : not used.
  
  // Function signature needs to be conformant with what integrate_ode_* accepts.
  real[] dSLD_dt(real t, real[] sld, real[] theta, real[] x_r, int[] x_i) {
    real dsld_dt = theta[1]*sld[1] - theta[2] *exp(-theta[3]*t)*theta[4]*log(1+x_r[1])* sld[1];
    return { dsld_dt};
  }
}

data {
  int<lower = 1> T;// total number of observations for all subjects ( without the initial time point (baseline) )
  int<lower=1> J;// total number of subjects
  int<lower=0> Nobs[J];// number of observations for each subjects without the baseline
  
  real ts[sum(Nobs)];// all unique time points  for which ODE solutions are requested
  real sld[T, 1]; // response vector of dimemnsion: (T,1); length of rows is the total number of observations T and 1 column
  real dose[T]; // dose amount at each observation
}
transformed data{ 
  real t0; // initial time point 
  real x_r[1]; // real data which in this case is the continous  covariate DOSE 
  real sigma;  // sd error 
  int x_i[0]; // ineteger covariate that we do not use here, it is a a vector of length 0
  
  
  int f_index[J];// create an index of leght total number of subjects J. We need this index later 
  f_index[1] =1;// initialize the thendex to 1 for the first subject
  for (i in 2:J)
    f_index[i] = f_index[i-1] + Nobs[i-1];
  x_r[1]=dose[T];// conitinous covariate dos e
  sigma = 0.17;
}
parameters {
  
  // define the parameters
  // each parameters should be estimated for each subjects
  // in particular we should have in teh end an estimate of the parameter for each subjects
  real<lower=0>kg[J];
  real<lower=0>ks0[J];
  real<lower=0>gamma[J];
  real<lower=0>alpha[J];
  real<lower=0>base[J,1];
  
}


 transformed parameters {
   
}


model{
  
  
  
  for (j in 1:J){
 base[j,1] ~ normal(60,0.1);
 }
kg ~ normal(0.15, 0.1);
ks0 ~ normal(0.15, 0.1);
gamma ~ normal(0.15, 0.1);
alpha ~ normal(1.5, 0.1);

}

generated quantities{
real theta[J,4];
real IPRED[sum(Nobs)];
  real DV[sum(Nobs)]; 
  

for(j in 1:J){ 
  
theta[j,1]=kg[j];
theta[j,2]=ks0[j];
theta[j,3]=gamma[j];
theta[j,4]=alpha[j]; 
}



for (i in 1:J) {
    real sol[Nobs[i],1];
   
 sol = integrate_ode_bdf(dSLD_dt, base[i] , 0.0, segment(ts, f_index[i], Nobs[i]), theta[i], x_r, x_i); 
 for (j in 1:Nobs[i])
     IPRED[f_index[i] + j-1] =sol[j,1];
}

for (j in 1:sum(Nobs))
// likelihood 
    DV[j] = IPRED[j]+normal_rng(0,sigma);
 
  
}

