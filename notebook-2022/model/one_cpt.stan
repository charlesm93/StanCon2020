
data {
  int N;
  vector[N] y_obs;
  vector[N] time;
  int n_patients;
  int<lower = 1, upper = N> start[n_patients];
  int<lower = 1, upper = N> end[n_patients];
  
  real<lower = 0> y0;  // initial dose.
}

parameters {
  real<lower = 0> sigma;
  real<lower = 0> sigma_0;
  real<lower = 0> sigma_1;
  real k_0_pop;
  real k_1_pop;
  vector[n_patients] k_0;
  vector[n_patients] k_1;
}

transformed parameters {
  vector[N] y;
  for (i in 1:n_patients) {
    for (j in start[i]:end[i]) {
      y[j] = y0 / (k_0[i] - k_1[i]) * k_1[i] 
        * (exp(- k_1[i] * time[j]) - exp(- k_0[i] * time[j]));
    }
  }
}

model {
  // priors
  sigma ~ normal(0, 1);
  sigma_0 ~ normal(0, 1);
  sigma_1 ~ normal(0, 1);
  k_0_pop ~ normal(2, 0.5);
  k_1_pop ~ normal(1, 0.5);
  k_0 ~ normal(k_0_pop, sigma_0);
  k_1 ~ normal(k_1_pop, sigma_1);

  // likelihood
  y_obs ~ normal(y, sigma);
}

generated quantities {
  real y_pred[N] = normal_rng(y, sigma);
}
