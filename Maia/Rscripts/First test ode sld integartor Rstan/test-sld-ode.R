library(rstan)

expose_stan_functions("sld-ode.stan")

dsld_dt(t = 5, y = 3, theta = c(0.1, 0.2, 0.3, 0.4), x_r = numeric(), x_i = integer())
dsld_dt(t = 5, y = 4, theta = c(0.1, 0.2, 0.3, 0.4), x_r = numeric(), x_i = integer())

sld(t = c(1, 2, 3), theta = c(0.1, 0.2, 0.3, 0.4), y0 = 1, t0 = 0, x_r = numeric(), x_i = integer())
