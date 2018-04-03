// Based on code from Michael Betancourt
// Defines a function for sampling from the predictive distribution, for visualization.
functions {
  vector gp_pred_rng(real[] x2_1, real[]x2_2,
                     vector y1, real[] x1_1,real[] x1_2,
                     real alpha, real rho_1,real rho_2, real sigma, real delta) {
    int N1 = rows(y1);
    int N2 = size(x2_1);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   cov_exp_quad(x1_1, alpha, rho_1).*cov_exp_quad(x1_2,1.,rho_2)
                         + diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_exp_quad(x1_1, x2_1, alpha, rho_1).*cov_exp_quad(x1_2,x2_2,1.,rho_2);
      vector[N2] f2_mu = (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 =   cov_exp_quad(x2_1, alpha, rho_1).*cov_exp_quad(x2_2,1.,rho_2) - v_pred' * v_pred
                              + diag_matrix(rep_vector(delta, N2));
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  //Data
  int<lower=1> N;
  real x_1[N];
  real x_2[N];
  vector[N] y;

  //Locations for predictions
  int<lower=1> N_predict;
  real x_1_predict[N_predict];
  real x_2_predict[N_predict];
}

parameters {
  real<lower=0> rho_1;
  real<lower=0> rho_2;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  // sets up squared exponential kernel, plus gaussian noise.
  matrix[N, N] cov =   cov_exp_quad(x_1, alpha, rho_1).*cov_exp_quad(x_2,1.,rho_2)
                     + diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  // priors on hyperparameters -- feel free to change!
  rho_1 ~ inv_gamma(1,1);
  rho_2 ~ inv_gamma(1,1);
  alpha ~ inv_gamma(1,1);
  sigma ~ inv_gamma(1,1);
  //likelihood
  y ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
}

generated quantities {
  // generates a random function
  vector[N_predict] f_predict = gp_pred_rng(x_1_predict, x_2_predict, y, x_1, x_2, alpha, rho_1, rho_2, sigma, 1e-10);

  //adds on noise
  vector[N_predict] y_predict;
  for (n in 1:N_predict)
    y_predict[n] = normal_rng(f_predict[n], sigma);
}

