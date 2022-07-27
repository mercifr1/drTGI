// TGI model as an hierarchical ODE : for 3 subjects
//  Estimating System Parameters and Initial State with rstan


functions {
  real[] dSLD_dt(real t, real[] sld, real[] theta, real[] x_r, int[] x_i) {
    real dsld_dt = theta[1]*sld[1] - theta[2] *exp(-theta[3]*t)*theta[4]*log(1+x_r[1])* sld[1];
    return { dsld_dt};
  }
}

transformed data {
  int<lower = 1> T=18;//number of observations 
  int<lower=1> J=3;// number of patients
  int<lower=0> Nobs[J]={7,7,4};// number of observations for each subjects
  real t0=0;
  real ts[sum(Nobs)]={0.122, 0.235, 0.360 ,0.467, 0.584, 0.691, 0.812,
  0.113, 0.223, 0.342, 0.441, 0.563,0.704, 0.795,
  0.124, 0.243, 0.345, 0.456};// time points to solve the ODE
  real sld[T, 1]={{68.60061}, {71.54244}, {75.75315}, 
  {79.42463}, {84.41742}, {88.65481}, {94.14883},
  {65.20813},{ 66.87848}, {69.10894}, {71.06266}, 
  {74.04987} ,{77.75904}, {80.49416},{76.62179}, {79.19184}, {81.60073}, {85.23003}};// response values sld 
  real dose[T]={1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67,
  1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67,5,5,5,5};// real covariate dose 
  real x_r[1];
  real sigma;
  int x_i[0];
  int f_index[J];
  f_index[1] =1;
  for (i in 2:J)
    f_index[i] = f_index[i-1] + Nobs[i-1];
  x_r[1]=dose[T];
  sigma=0.17;
}
parameters {
  real<lower=0> kg[J];
  real<lower=0> ks0[J];
  real<lower=0> gamma[J];
  real<lower=0> alpha[J];
  real<lower=0> base[J];
  
  
  
  


}

transformed parameters{
  real theta[4];
  real y0[1];
real IPRED[sum(Nobs)];
for(j in 1:J){ 
  y0[1]=exp(log(base[j]));
theta[1]=kg[j];
theta[2]=ks0[j];
theta[3]=gamma[j];
theta[4]=alpha[j]; 
}



for (i in 1:J) {
    real sol[Nobs[i],1];
   
 sol = integrate_ode_bdf(dSLD_dt, y0 , 0.0, segment(ts, f_index[i], Nobs[i]), theta, dose, x_i); 
 for (j in 1:Nobs[i])
     IPRED[f_index[i] + j-1] =sol[j,1];
}

}



model{
  
    

 
 
base ~ normal(80,10);
kg ~ normal(0.5, 0.5);
ks0 ~ normal(0.5, 0.5);
gamma ~ normal(0.5, 0.5);
alpha ~ normal(0, 2);


// likelihood 

    sld[,1] ~ normal(IPRED,sigma);
}

generated quantities{


 
  
}



