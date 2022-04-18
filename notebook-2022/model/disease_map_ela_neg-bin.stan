
functions {
  matrix K_functor (vector[] x, int n_obs, real alpha, real rho) {
    matrix[n_obs, n_obs] K = cov_exp_quad(x, alpha, rho);
    for (i in 1:n_obs) K[i, i] += 1e-8;
    return K;
  }
  
  // this must follow a strict signature
  real L_functor(
    vector theta,        // latent Gaussian variable
    vector eta,          // hyperparameters for the likelihood
    vector log_ye,       // real data  (mean offset)
    int[] y) {           // interger data (observations)
      int n = num_elements(theta);
      real alpha[n] = to_array_1d(log_ye + theta);
      vector[n] eta_rep = rep_vector(eta[1], n);
      return neg_binomial_2_lpmf(y | exp(alpha), eta_rep);
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
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  
  // arguments for likelihood functor
  vector[n_obs] log_ye = log(ye);
  int delta_int[n_obs + 1];
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
  real<lower = 0> eta;
}

transformed parameters {
  vector[1] eta_vec = to_vector({eta});
}

model {
  rho ~ inv_gamma(rho_location_prior, rho_scale_prior);
  alpha ~ inv_gamma(alpha_location_prior, alpha_scale_prior);
  eta ~ normal(0, 1);

  target += laplace_marginal_lpmf(y | L_functor, eta_vec, log_ye,
                                      theta_0,
                                      K_functor, x, n_obs, alpha, rho);
}

generated quantities {
  vector[n_obs] theta = laplace_marginal_rng(L_functor, eta_vec, log_ye, y,
                                             theta_0, K_functor,
                                             forward_as_tuple(x, n_obs),
                                             forward_as_tuple(x, n_obs),
                                             alpha, rho);
}
