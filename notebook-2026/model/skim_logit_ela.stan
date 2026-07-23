// same as skim_pairwise_logit_ela.stan model, but rewritten
// to be compatible with ovarian and prostate data set.

functions {

  real likelihood_function (vector theta, array[] int y, int x_dummy) {
    return bernoulli_logit_lpmf(y | theta);
  }

  matrix covariance_function (vector kappa_squared,
                              real eta_one,
                              real eta_two,
                              matrix x,
                              matrix x2,
                              int n,
                              int d,
                              real scale_icept) {
    matrix[n, n] K1 = diag_post_multiply(x, kappa_squared) * x';
    matrix[n, n] K2 = diag_post_multiply(x2, kappa_squared) * x2';
    matrix[n, n] K = 0.5 * square(eta_two) * square(K1 + 1) -
      0.5 * square(eta_two) * K2 +
      (square(eta_one) - square(eta_two)) * K1 +
      square(scale_icept) - 0.5 * square(eta_two);

    // add jitter
    return add_diag(K, 1e-5);;
  }
}

data {
  int<lower = 1> n;
  int<lower = 1> d;
  array[n] int y;
  matrix[n, d] x;

  // variables for priors
  real<lower = 0> scale_icept;
  real<lower = 0> scale_global;
  real<lower = 1> nu_global;
  real<lower = 1> nu_local;
  real<lower = 0> slab_scale;
  real slab_df;
}

transformed data {
   real d0 = 5;  // expected number of active predictors

  // variables for prior
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5 * slab_df;
  real sigma = 4;

  vector[n] mu = rep_vector(0, n);
  matrix[n, d] x2 = square(x);

  // variables for built-in likelihood
  // array[n] int y_index;
  // for (i in 1:n) y_index[i] = i-1;
  // vector[n] m = rep_vector(0, n);
  int hessian_block_size = 1;
}

parameters {
  real<lower = 0> tau;
  vector<lower = 0>[d] lambda;
  real<lower = 0> caux;
  real<lower = 0> xi;
}

transformed parameters {
  real<lower = 0> eta_one = (d0 / (d - d0))
    * (sigma / sqrt(n)) * tau;
  real<lower = 0> m_squared = slab_scale2 * caux;
  vector[d] kappa_squared = m_squared * square(lambda)
    ./ (m_squared + square(eta_one) * square(lambda));
  real<lower = 0> eta_two = square(eta_one) / m_squared * xi;

  vector[d + 2] parm;
  parm[1:d] = lambda;
  parm[d + 1] = eta_two;
  parm[d + 2] = eta_one;
}

model {
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, 1);
  caux ~ inv_gamma(half_slab_df, half_slab_df);
  xi ~ inv_gamma(half_slab_df, half_slab_df);

  // using user-defined function
  target += laplace_marginal(likelihood_function, (y, 0), hessian_block_size,
        covariance_function,
        (kappa_squared, eta_one, eta_two, x, x2, n, d, scale_icept));

  // using built-in function
  // y ~ laplace_marginal_bernoulli_logit(y_index, m, hessian_block_size,
  //       covariance_function,
  //       (kappa_squared, eta_one, eta_two, x, x2, n, d, scale_icept));
}

generated quantities {
  vector[n] theta = laplace_latent_rng(likelihood_function, (y, 0),
        hessian_block_size, covariance_function,
        (kappa_squared, eta_one, eta_two, x, x2, n, d, scale_icept));

  // vector[n] theta = laplace_latent_bernoulli_logit_rng(y, y_index, m,
  //       hessian_block_size, covariance_function,
  //       (kappa_squared, eta_one, eta_two, x, x2, n, d, scale_icept));

  vector[n] p = inv_logit(theta);
}


