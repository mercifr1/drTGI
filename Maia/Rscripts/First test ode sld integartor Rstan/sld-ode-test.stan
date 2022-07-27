functions{
  // t: single time point at which we want to evaluate the derivative.
  // y: the scalar function value.
  // theta: parameter vector (4 parameters).
  // x_r, x_i: not used.
  // Function signature needs to be conformant with what integrate_ode_* accepts.
  real[] dsld_dt(real t, real[] y, real[] theta, real[] x_r, int[] x_i){

    real dydt[1];

    // Define the derivative that defines the SLD.
    dydt[1] = theta[1] * y[1] - theta[2] * exp(-theta[3] * t) * theta[4] *log(1+x_r[1])* y[1];
    return dydt;
  }

  // t: multiple time points at which we want to evaluate the SLD function.
  // theta: parameter vector (4 parameters).
  // y0: function value at time t0.
  // t0: initial time point (baseline).
  real[] sld(real[] t, real[] theta, real y0, real t0, data real[] x_r, data int[] x_i) {
    real y[1] = {y0};

    // see https://mc-stan.org/docs/functions-reference/functions-old-ode-solver.html
    real result[size(t), 1] = integrate_ode_rk45(dsld_dt, y, t0, t, theta, x_r, x_i);
    return result[:, 1];
  }
}

data {
  real x_r[1];
  int x_i[0];
}









