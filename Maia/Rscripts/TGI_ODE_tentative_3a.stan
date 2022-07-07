functions {
real[] sld(real t, // time
real[] y, // state
real[] theta, // parameters
real[] x_r,  // data (real)
int[] x_i) // data(integer)
{
real dydt[1];
dydt[1] = theta[1] * y[1]-theta[2]*exp(-t*theta[3])*theta[4]*log(1+x_r[1])*y[1];
return dydt;
  }
}


data {
int<lower=1> N; // total nr. of observations
int<lower=1> J; // nr. of subjects
int<lower=1> T; // nr. of timepoints (all unique timepoints)

real t0; // starting time
real<lower=t0> ts[T]; // timepoints for which ODE solutions are requested

real y[N,1]; // response vector
 real dose; // dose amount
int<lower=1> nTime[J]; // nr. of observations per subject
int<lower=1> timeInd[N]; // numeric time indicator

}
transformed data {
real x_r[1];
int x_i[0];
x_r[1]=dose;
}
parameters {
  real kg[J];
  real ks0[J];
  real gamma[J];
  real alpha[J];
  real base[J];
  real<lower=0> sigma;
  
  
real<lower=0> mu_kg;
real<lower=0> sigma_kg_sq;
real<lower=0> mu_ks0;
real<lower=0> sigma_ks0_sq;
real<lower=0> mu_gamma;
real<lower=0> sigma_gamma_sq;
real<lower=0> mu_alpha;
real<lower=0> sigma_alpha_sq;
real<lower=0> mu_base;
real<lower=0> sigma_base_sq;

}
transformed parameters {

}
model {
real y_hat[T,1];
int index;
//prior distributions
sigma ~ normal(0,5);
mu_kg ~ lognormal(3.61, 0.002);
sigma_kg_sq ~ inv_gamma(54.441, 24.476);

mu_ks0 ~ lognormal(-0.115, 0.514);
sigma_ks0_sq ~ inv_gamma(54.441, 24.476);

mu_gamma ~ lognormal(0.2,0.05);
sigma_gamma_sq ~ inv_gamma(50,20);

mu_alpha ~ lognormal(0.2,0.05);
sigma_alpha_sq ~ inv_gamma(50,20);

mu_base ~ lognormal(0.2,0.05);
sigma_base_sq ~ inv_gamma(50,20);


index = 0;



for(j in 1:J){
real y0[1];
real theta[4];
//initial conditions
 y0[1] = base[j];

theta[1] = kg[j];
theta[2] = ks0[j];
theta[3] = gamma[j];
theta[4] = alpha[j];

y_hat = integrate_ode_rk45(sld, y0, t0, ts, theta, x_r, x_i);
//likelihood
for (t in 1:nTime[J]) {
index = index + 1;
y[index] ~ normal(y_hat[timeInd[index]], sigma);}
base[j] ~ normal(mu_base, sigma_base_sq);
kg[j] ~ normal(mu_kg, sigma_kg_sq);
ks0[j] ~ normal(mu_ks0, sigma_ks0_sq);
gamma[j] ~ normal(mu_gamma, sigma_gamma_sq);
alpha[j] ~ normal(mu_alpha, sigma_alpha_sq);}
}