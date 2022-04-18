// same as skim_pairwise_logit.stan model, but rewritten
// to be compatible with ovarian and prostate data set.

data {
  int<lower = 0> n;
  int<lower = 0> d;
  int<lower = 0, upper = 1> y[n];
  matrix[n, d] x;

  // variables for priors
  real<lower = 0> scale_icept;
  real<lower = 0> scale_global;
  real<lower = 1> nu_global;
  real<lower = 1> nu_local;
  real<lower = 0> slab_scale;
  real<lower = 0> slab_df;
}

transformed data {
   real d0 = 5;  // expected number of active predictors

  // variables for prior
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5 * slab_df;

  vector[n] mu = rep_vector(0, n);
  matrix[n, d] x2 = square(x);
}

parameters {
  vector[n] z;
  real<lower = 0> tau;
  vector<lower = 0>[d] lambda;
  real<lower = 0> caux;
  real<lower = 0> xi;
}

transformed parameters {
  real<lower = 0> eta_one = scale_global * tau;
  real<lower = 0> m_squared = slab_scale2 * caux;
  vector[d] kappa_squared = m_squared * square(lambda)
    ./ (m_squared + square(eta_one) * square(lambda));
  real<lower = 0> eta_two = square(eta_one) / m_squared * xi;

  vector[n] f;
  {
    matrix[n, n] L_K;
    matrix[n, n] K1 = diag_post_multiply(x, kappa_squared) * x';
    matrix[n, n] K2 = diag_post_multiply(x2, kappa_squared) * x2';
    matrix[n, n] K = 0.5 * square(eta_two) * square(K1 + 1) -
      0.5 * square(eta_two) * K2 +
      (square(eta_one) - square(eta_two)) * K1 +
      square(scale_icept) - 0.5 * square(eta_two);

    // add jitter term for numerical stability
    for (i in 1:n) K[i, i] += 1e-5;
    L_K = cholesky_decompose(K);
    f = mu + L_K * z;
  }
}

model {
  // halft-t priors for lambdas and tau, and inverse-gamma for c^2
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, 1);
  caux ~ inv_gamma(half_slab_df, half_slab_df);
  xi ~ inv_gamma(half_slab_df, half_slab_df);

  z ~ normal(0, 1);
  y ~ bernoulli_logit(f);
}

generated quantities {
  vector[n] p = inv_logit(f);
}
