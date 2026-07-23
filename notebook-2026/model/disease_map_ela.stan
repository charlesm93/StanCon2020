
functions {
  real likelihood_function (vector theta, vector log_ye, array[] int y) {
    int n = num_elements(theta);
    vector[n] alpha = log_ye + theta;
    return sum(to_vector(y) .* alpha - exp(alpha));
  }

  matrix cov_function (array[] vector x, int n_obs, real alpha, real rho) {
    matrix[n_obs, n_obs] K = gp_exp_quad_cov(x, alpha, rho);
    return add_diag(K, 1e-8);
  }
}

data {
  int n_obs;
  int n_coordinates;
  array[n_obs] int y;
  vector[n_obs] ye;
  array[n_obs] vector[n_coordinates] x;
}

transformed data {
  vector[n_obs] log_ye = log(ye);
  int hessian_block_size = 1;
  
  real tol = 1e-6;
  int max_num_steps = 100;
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  
  // variables for built-in likelihood
  // array[n_obs] int y_index = linspaced_int_array(n_obs, 1, n_obs);
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
}

model {
  rho ~ inv_gamma(2.4, 14.8);
  alpha ~ inv_gamma(10, 10);

  target += laplace_marginal(likelihood_function, (log_ye, y),
                             hessian_block_size,
                             cov_function, (x, n_obs, alpha, rho));
  
  // alternative using built-in function
  // y ~ laplace_marginal_poisson_log(y_index, log_ye, hessian_block_size,
  //                                  cov_function, (x, n_obs, alpha, rho));
}

generated quantities {
  vector[n_obs] theta
    = laplace_latent_rng(likelihood_function, (log_ye, y),
                         hessian_block_size,
                         cov_function, (x, n_obs, alpha, rho));

  // alternative using built-in function
  // vector[n_obs] theta
  //   = laplace_latent_poisson_log_rng(y, y_index, log_ye, hessian_block_size,
  //                                    cov_function, (x, n_obs, alpha, rho));
}
