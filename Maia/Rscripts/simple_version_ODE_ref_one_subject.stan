functions {
  real[] dSLD_dt(real t, real[] sld, real[] theta, real[] x_r, int[] x_i) {
    real dsld_dt = theta[1]*sld[1] - theta[2] *exp(-theta[3]*t)*theta[4]*log(1+x_r[1])* sld[1];
    return { dsld_dt};
  }
}

transformed data {
  int<lower = 1> T=7;
  real t0=0;
  real ts[T]={0.122,0.235,0.360,0.467,0.584,0.691,0.812};
  real sld[T, 1]={{68.60061}, {71.54244}, {75.75315}, {79.42463}, {84.41742}, {88.65481}, {94.14883}};
  real dose[T]={1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67};
  real x_r[1];
  
  int x_i[0];
  x_r[1]=dose[T];
}
parameters {
  real<lower=0> kg;
  real<lower=0> ks0;
  real<lower=0> gamma;
  real<lower=0> alpha;
  real<lower=0> y0[1];
  real sigma;
  
  


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



model {
  
  real y_hat[T,1];
  real theta[4];
sigma ~ normal(0.1, 0.1);
y0 ~ normal(60, 0.1);
kg ~ normal(0.15, 0.1);
ks0 ~ normal(0.15, 0.1);
gamma ~ normal(0.15, 0.1);
alpha ~ normal(1.5, 0.1);
theta[1]=kg;
theta[2]=ks0;
theta[3]=gamma;
theta[4]=alpha;
 y_hat = integrate_ode_bdf(dSLD_dt, y0 , 0.0, ts, theta, x_r, x_i); 
 
 sld[,1] ~ normal(y_hat[, 1],  sigma);
  
}
generated quantities {
}

