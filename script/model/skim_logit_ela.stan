// same as skim_pairwise_logit_ela.stan model, but rewritten
// to be compatible with ovarian and prostate data set.

functions {
  matrix K_functor (vector parm, matrix x_tot, 
                    real[] delta, int[] delta_int) {
    int n = delta_int[1];
    int d = delta_int[2];
    real scale_icept = delta[1];
    vector[d] kappa_squared = parm[1:d];
    real eta_two = parm[d + 1];
    real eta_one = parm[d + 2];

    matrix[n, n] x;
    matrix[n, n] x2;
    matrix[n, n] K1;
    matrix[n, n] K2;
    matrix[n, n] K;

    x = x_tot[1:n, ];
    x2 = x_tot[(n + 1):(2 * n), ];
    K1 = diag_post_multiply(x, kappa_squared) * x';
    K2 = diag_post_multiply(x2, kappa_squared) * x2';
    K = 0.5 * square(eta_two) * square(K1 + 1) -
      0.5 * square(eta_two) * K2 +
      (square(eta_one) - square(eta_two)) * K1 +
      square(scale_icept) - 0.5 * square(eta_two);

    // add jitter
    for (i in 1:n) K[i, i] += 1e-5;

    return K;
  }
}

data {
  int<lower = 1> n;
  int<lower = 1> d;
  int y[n];
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
  
  // variables for the Laplace approximation
  real delta[1] = {scale_icept};
  int delta_int[2] = {n, d};
  int n_samples[n] = rep_array(1, n);
  vector[n] theta0 = rep_vector(0, n);
  matrix[2 * n, d] x_tot;
  x_tot[1:n, ] = x;
  x_tot[(n + 1):(2 * n), ] = square(x);
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
  
  target += 
    laplace_marginal_bernoulli_logit_lpmf(y | n_samples, K_functor, parm,
                                          x_tot, delta, delta_int, theta0);
}

generated quantities {
  vector[n] p = inv_logit(
    laplace_bernoulli_logit_rng(y, n_samples, K_functor, parm,
                                x_tot, delta, delta_int, theta0));
}
