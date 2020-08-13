
functions {
  matrix K (vector phi, vector[] x, real[] delta, int[] delta_int) {
    real alpha = phi[1];
    real rho = phi[2];
    return add_diag(cov_exp_quad(x, alpha, rho), 1e-8);
  }
}

data {
  int n_obs;
  int n_coordinates;
  int y[n_obs];
  vector[n_obs] ye;
  vector[n_coordinates] x[n_obs];

  real<lower = 0> alpha_location_prior;
  real<lower = 0> alpha_scale_prior;
  real<lower = 0> rho_location_prior;
  real<lower = 0> rho_scale_prior;
}

transformed data {
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  real delta[0];
  int delta_int[0];
  int n_samples[n_obs] = rep_array(1, n_obs);
  int n_phi = 2;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
}

transformed parameters {
  vector[n_phi] phi = to_vector({alpha, rho});
}

model {
  alpha ~ inv_gamma(alpha_location_prior, alpha_scale_prior);
  rho ~ inv_gamma(rho_location_prior, rho_scale_prior);

  // y ~ laplace_marginal_poisson_log(n_samples, ye, K, phi, x, delta_int, theta_0);

  target += laplace_marginal_poisson_log_lpmf(y | n_samples, ye, K, phi,
                                         x, delta, delta_int, theta_0);
}

generated quantities {
  vector[n_obs] theta = laplace_poisson_log_rng(y, n_samples, ye, K, phi, x,
                                                delta, delta_int, theta_0);
}
