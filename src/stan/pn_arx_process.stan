//
// This Stan program fits the projected normal AR(1) model

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, 2] U;
  matrix[N, K] X;
  real<lower=0> sigma_0;
  
  int<lower=0> N_ppd;
  matrix[N_ppd, K] X_ppd;
  
  int<lower=0, upper=1> weakly_inf_prior;
}

transformed data {
  vector[K] X_bar = 1.0 / N * (rep_vector(1.0, N)' * X)';
  matrix[N, K] X_centered = X - rep_matrix(X_bar', N);
  X_centered[,1] = rep_vector(1, N);
  
  vector[K] X_ppd_bar = 1.0 / N_ppd * (rep_vector(1.0, N_ppd)' * X_ppd)';
  matrix[N_ppd, K] X_ppd_centered = X_ppd - rep_matrix(X_ppd_bar', N_ppd);
  X_ppd_centered[,1] = rep_vector(1, N_ppd);
}

parameters {
  vector[2] mu_0;
  matrix[K,2] B_mat;
  matrix<lower=-1,upper=1>[2,2] auto_cor_mat;
  vector<lower=0>[N] latent_lengths;
  
  real<lower=0> sigma_w;
  real<lower=-1, upper=1> rho_w;
}

transformed parameters {
  matrix[N, 2] Y;
  
  Y = diag_pre_multiply(latent_lengths, U);
  
  matrix[2,2] Sigma_inv;
  Sigma_inv[1,1] = 1;
  Sigma_inv[1,2] = -sigma_w * rho_w;
  Sigma_inv[2,1] = -sigma_w * rho_w;
  Sigma_inv[2,2] = 1.0 / sigma_w^2;
  Sigma_inv = 1 / (sigma_w^2 * (1 - rho_w^2)) * Sigma_inv;
}

model {
  matrix[N,2] mu;
  latent_lengths ~ normal(0, 10);
  sigma_w ~ inv_chi_square(3);
  rho_w ~ normal(0, 1);
  
  if (weakly_inf_prior == 1) {
    mu_0 ~ normal(0, 10);
    to_vector(B_mat) ~ normal(0, 10);
    to_vector(auto_cor_mat) ~ normal(0, 1);
  }
  
  mu[1,] = (mu_0' + X_centered[1,] * B_mat);
  real A = quad_form(Sigma_inv, U[1,]');
  real B = U[1,] * Sigma_inv * mu[1,]';
  target += log(latent_lengths[1]) - 1.0 / 2 * A * (latent_lengths[1] - B / A)^2;
  
  Y[1,] ~ multi_normal(mu[1,], Sigma_inv);
  
  // Latent Length sampling
  for (t in 2:N) {
      mu[t,] = (mu_0' + Y[t-1,] * auto_cor_mat + X_centered[t,] * B_mat);
      A = quad_form(Sigma_inv, U[t,]');
      B = U[t,] * Sigma_inv * mu[t,]';
      target += log(latent_lengths[t]) - 1.0 / 2 * A * (latent_lengths[t] - B / A)^2;
      Y[t,] ~ multi_normal_prec(mu[t,], Sigma_inv);
      
  }
  
  // // Sampling for Y
  // for (t in 1:N) {
  // }
}

generated quantities {// Post-processed Parameters
  
  matrix[2,2] Sigma;
  Sigma[1,1] = sigma_w^2;
  Sigma[1,2] = sigma_w * rho_w;
  Sigma[2,1] = sigma_w * rho_w;
  Sigma[2,2] = 1;
  // to get the constrained version
  // where sigma_w2i = 1
  // vector<lower=0>[2] c_diag = rep_vector(Sigma[2, 2]^-0.5, 2);
  // 
  // matrix[2, 2] Sigma_w = quad_form_diag(Sigma, c_diag);
  // matrix[K, 2] B_c = B_mat * diag_matrix(c_diag);
  // matrix[2, 2] auto_cor_mat_c = diag_post_multiply(auto_cor_mat, c_diag);
  // vector[2] mu_0_c = diag_matrix(c_diag) * mu_0;
  // 
  // // Variance Components for Angular Components
  // real<lower=0> sigma_w = Sigma_w[1,1]^0.5;
  // real<lower=-1, upper=1> rho_w = Sigma_w[1,2] / sigma_w;
  
  real<lower=-pi(), upper=pi()> theta = atan2(mu_0[2], mu_0[1]);
  
  // Draws for (1,...,T,T+1,...,T+T_ppd) (Not one-shot prediction)
  matrix[N+N_ppd, K] X_mat = append_row(X_centered, X_ppd_centered);
  matrix[N+N_ppd, 2] mu_ppd;
  mu_ppd[1,] = mu_0';
  
  matrix[N+N_ppd,2] eta = X_mat * B_mat;
  
  matrix[N+N_ppd, 2] Y_ppd;
  Y_ppd[1,] = multi_normal_rng(mu_ppd[1,] + eta[1,], Sigma)';
  
  // ppd draws projected to S^1
  matrix[N+N_ppd, 2] U_ppd;
  U_ppd[1,] = Y_ppd[1,] / sqrt(dot_self(Y_ppd[1,]));
  
   // ppd draws converted to angles
   vector<lower=-pi(), upper=pi()>[N+N_ppd] theta_ppd;
   theta_ppd[1] = atan2(U_ppd[1,2], U_ppd[1,1]);
   
  for (t in 2:(N+1)) {
    mu_ppd[t,] = (mu_0 + auto_cor_mat * Y[t-1,]')';
    // Prediction/backcast for time t given observed Y_1,...,Y_{N+1}
    Y_ppd[t, ] = multi_normal_rng(mu_ppd[t,] + eta[t,], Sigma)';
    U_ppd[t, ] = Y_ppd[t,] / sqrt(dot_self(Y_ppd[t,]));
    theta_ppd[t] = atan2(U_ppd[t, 2], U_ppd[t, 1]);
  }
  
  for (t in 2:N_ppd) {
    mu_ppd[N+t,] = (mu_0 + auto_cor_mat * Y_ppd[N+t-1,]')';
    // Forecast for time t given observed Y_N+2,...,Y_{N+N_ppd}
    Y_ppd[N+t, ] = multi_normal_rng(mu_ppd[N+t,] + eta[N+t,], Sigma)';
    U_ppd[N+t, ] = Y_ppd[N+t,] / sqrt(dot_self(Y_ppd[N+t,]));
    theta_ppd[N+t] = atan2(U_ppd[N+t, 2], U_ppd[N+t, 1]);
  }
  
}
