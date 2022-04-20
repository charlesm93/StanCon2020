
functions {
  matrix K_functor (vector[] x, int n_obs, real alpha, real rho) {
    matrix[n_obs, n_obs] K = cov_exp_quad(x, alpha, rho);
    for (i in 1:n_obs) K[i, i] += 1e-8;
    return K;
  }
}

data {
  int n_obs;
  int n_coordinates;
  int y[n_obs];
  vector[n_obs] ye;
  array[n_obs] vector[n_coordinates] x;
  real rho_location_prior;
  real rho_scale_prior;
  real alpha_location_prior;
  real alpha_scale_prior;
}

transformed data {
  real tol = 1e-6;
  int max_num_steps = 100;
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  int n_samples[n_obs] = rep_array(1, n_obs);
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
}

model {
  rho ~ inv_gamma(rho_location_prior, rho_scale_prior);
  alpha ~ inv_gamma(alpha_location_prior, alpha_scale_prior);

  target += laplace_marginal_poisson_2_log_lpmf(y | n_samples, ye, theta_0, K_functor,
                                                x, n_obs, alpha, rho);
}

generated quantities {
  vector[n_obs] theta
    = laplace_marginal_poisson_2_log_rng(y, n_samples, ye, theta_0, K_functor,
                                         forward_as_tuple(x, n_obs), 
                                         forward_as_tuple(x, n_obs), 
                                         alpha, rho);
}
