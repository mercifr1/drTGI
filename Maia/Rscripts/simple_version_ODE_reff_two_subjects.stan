functions {
  real[] dSLD_dt(real t, real[] sld, real[] theta, real[] x_r, int[] x_i) {
    real dsld_dt = theta[1]*sld[1] - theta[2] *exp(-theta[3]*t)*theta[4]*log(1+x_r[1])* sld[1];
    return { dsld_dt};
  }
}

transformed data {
  int<lower = 1> T=14;//number of observations
  int<lower=1> J=2;// number of patients
  int<lower=0> Nobs[J]={7,7};// number of observations for each subjects
  real t0=0;
  real ts[sum(Nobs)]={0.122, 0.235, 0.360 ,0.467, 0.584, 0.691, 0.812,  0.113, 0.223, 0.342, 0.441, 0.563,0.704, 0.795};
  real sld[T, 1]={{68.60061}, {71.54244}, {75.75315}, 
  {79.42463}, {84.41742}, {88.65481}, {94.14883},
  {65.20813},{ 66.87848}, {69.10894}, {71.06266}, 
  {74.04987} ,{77.75904}, {80.49416}};
  real dose[T]={1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67,
  1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67};
  real x_r[1];
  real sigma;
  int x_i[0];
  int f_index[J];
  f_index[1] =1;
  for (i in 2:J)
    f_index[i] = f_index[i-1] + Nobs[i-1];
  x_r[1]=dose[T];
  sigma = 0.17;
}
parameters {
  real<lower=0> kg;
  real<lower=0> ks0;
  real<lower=0> gamma;
  real<lower=0> alpha;
  real<lower=0> y0[1];
  //real sigma;
  
  


}


// transformed parameters {
//   
//   real<lower=0> theta[4];
//     real eta[5];
//     real<lower=0> y0[1];
//        y0[1]=exp(base+eta[1]*omega_base);
//     theta[1]=exp(kg+eta[2]*omega_kg);
//     theta[2]=exp(ks0+eta[3]*omega_ks0);
//     theta[3]= exp(gamma+eta[4]*omega_gamma);
//     theta[4]=exp(alpha+eta[5]*omega_alpha);
//     
//     
//    
// }


model{
  
  
 // sigma ~ normal(2, 0.2);
y0 ~ normal(60, 0.1);
kg ~ normal(0.15, 0.1);
ks0 ~ normal(0.15, 0.1);
gamma ~ normal(0.15, 0.1);
alpha ~ normal(1.5, 0.1);

}

generated quantities{
real theta[4];
real IPRED[sum(Nobs)];
  real DV[sum(Nobs)]; 
  
 theta[1]=kg;
theta[2]=ks0;
theta[3]=gamma;
theta[4]=alpha; 
  



for (i in 1:J) {
    real sol[Nobs[i],1];
 sol = integrate_ode_bdf(dSLD_dt, y0 , 0.0, segment(ts, f_index[i], Nobs[i]), theta, x_r, x_i); 
 for (j in 1:Nobs[i])
     IPRED[f_index[i] + j-1] =sol[j,1];
}

for (j in 1:sum(Nobs))
// likelihood 
    DV[j] = IPRED[j]+normal_rng(0,sigma);
 
  
}



