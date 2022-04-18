
functions {
    matrix K_functor(real sigma_0, real sigma_1, int n_patients) {
      real sigma_0_squared = sigma_0^2;
      real sigma_1_squared = sigma_1^2;
      vector[2 * n_patients] K_vec;
      int index_0;
      int index_1;

      for (i in 1:n_patients) {
        index_0 = 2 * (i - 1) + 1;
        index_1 = 2 * i;
        K_vec[index_0] = sigma_0_squared;
        K_vec[index_1] = sigma_1_squared;
      }

      return diag_matrix(K_vec);
    }

    real L_functor(vector theta, vector eta, vector delta, int[] delta_int) {
      int n_patients = delta_int[1];
      int N = delta_int[2];
      int start[n_patients] = delta_int[3:(n_patients + 2)];
      int end[n_patients] = delta_int[(n_patients + 3):(2 * n_patients + 2)];
      vector[N] y_obs = delta[1:N];
      vector[N] time = delta[(N + 1):(2 * N)];
      real y0 = delta[2 * N + 1];
      real k_0;
      real k_1;
      real k_0_pop = eta[1];
      real k_1_pop = eta[2];
      real sigma = eta[3];

      vector[N] y;

      for (i in 1:n_patients) {
        for (j in start[i]:end[i]) {
          k_0 = theta[2 * (i - 1) + 1] + k_0_pop;
          k_1 = theta[2 * i] + k_1_pop;

          y[j] = y0 / (k_0 - k_1) * k_1
            * (exp(- k_1 * time[j]) - exp(- k_0 * time[j]));
        }
      }

    return normal_lpdf(y_obs | y, sigma);
  }
}


data {
  int N;
  vector[N] y_obs;
  vector[N] time;
  int n_patients;
  array[n_patients] int<lower = 1, upper = N> start;
  array[n_patients] int<lower = 1, upper = N> end;
  real<lower = 0> y0;  // initial dose.
}

transformed data {
  // initial guess for Newton solver
  vector[2 * n_patients] theta0 = rep_vector(0, 2 * n_patients);

  // control parameters for Laplace marginal
  real tol = 1e-3;
  int max_num_steps = 1000;
  int hessian_block_size = 2;
  int solver = 3;
  int max_steps_line_search = 10;
  
  // arguments for L_functor
  array[2 * n_patients + 4] int delta_int;
  vector[2 * N + 1] delta;
  delta_int[1] = n_patients;
  delta_int[2] = N;
  delta_int[3:(n_patients + 2)] = start;
  delta_int[(n_patients + 3):(2 * n_patients + 2)] = end;
  delta[1:N] = y_obs;
  delta[(N + 1):(2 * N)] = time;
  delta[2 * N + 1] = y0;

}

parameters {
  real<lower = 0> sigma;
  real<lower = 0> sigma_0;
  real<lower = 0> sigma_1;
  real k_0_pop;
  real k_1_pop;
}

transformed parameters {
  vector[3] eta = to_vector({k_0_pop, k_1_pop, sigma});
}

model {
  // priors
  sigma ~ normal(0, 1);
  sigma_0 ~ normal(0, 1);
  sigma_1 ~ normal(0, 1);
  k_0_pop ~ normal(2, 0.5);
  k_1_pop ~ normal(1, 0.5);

  // likelihood
  target += laplace_marginal_tol_lpdf(delta | L_functor, eta, delta_int,
                                  tol, max_num_steps,
                                  hessian_block_size,
                                  solver, max_steps_line_search, theta0, K_functor,
                                  sigma_0, sigma_1, n_patients);
}

// generated quantities {
//   // vector[2 * n_patients] theta_pred = laplace_rng(L_f, eta, delta, delta_int,
//   //                                                 K_f, phi, x_mat_dummy,
//   //                                                 to_array_1d(delta), delta_int,
//   //                                                 theta0,
//   //                                                 tol, max_num_steps,
//   //                                                 hessian_block_size,
//   //                                                 solver);
//   //
//   // vector[n_patients] k_0;
//   // vector[n_patients] k_1;
//   //
//   // for (i in 1:n_patients) {
//   //   k_0[i] = k_0_pop + theta_pred[1 + 2 * (i - 1)];
//   //   k_1[i] = k_1_pop + theta_pred[2 * i];
//   // }
// }
