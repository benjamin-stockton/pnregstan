//
// This Stan program fits the projected normal AR(1) model

data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  array[N_obs] int ind_obs;
  array[N_mis] int ind_mis;
  int<lower=0> K;
  matrix[N_obs, 2] U_obs;
  matrix[N_mis+N_obs, K] X;
  real<lower=0> sigma_0;
  
  int<lower=0> N_ppd;
  matrix[N_ppd, K] X_ppd;
}

transformed data {
  int<lower=0> N = N_obs + N_mis;
  vector[K] X_bar = 1.0 / N * (rep_vector(1.0, N)' * X)';
  matrix[N, K] X_centered = X - rep_matrix(X_bar', N);
  
  // matrix[N, 2] U_lag = append_row(rep_vector(0,2), U[1:N,]);
  // X_centered = append_col(X_centered, U_lag);
  
  vector[K] X_ppd_bar = 1.0 / N_ppd * (rep_vector(1.0, N_ppd)' * X_ppd)';
  matrix[N_ppd, K] X_ppd_centered = X_ppd - rep_matrix(X_ppd_bar', N_ppd);
}

parameters {
  vector[2] mu_0;
  matrix[K,2] B_mat;
  matrix[2,2] auto_cor_mat;
  vector<lower=0>[N] latent_lengths;
  
  real<lower=0> sigma_w;
  real<lower=0, upper=1> pre_rho_w;
  
  matrix[N_mis, 2] Y_mis;
}

transformed parameters {
  real<lower=-1, upper=1> rho_w = 2*pre_rho_w - 1;
  matrix[N_obs, 2] Y_obs;
  matrix[N, 2] Y;
  
  Y_obs = diag_pre_multiply(latent_lengths[ind_obs], U_obs);
  
  matrix[2,2] Sigma;
  Sigma[1,1] = sigma_w^2;
  Sigma[1,2] = sigma_w * rho_w;
  Sigma[2,1] = sigma_w * rho_w;
  Sigma[2,2] = 1;
  
  Y[ind_obs,] = Y_obs;
  Y[ind_mis,] = Y_mis;
}

model {
  // mu_0 ~ normal(0, 100);
  // to_vector(auto_cor_mat) ~ normal(0, 100);
  // to_vector(B_mat) ~ normal(0, 100);
  latent_lengths ~ normal(0, 100);
  
  // sigma_w ~ inv_gamma(100, 0.01);
  // pre_rho_w ~ beta(1, 1);
  
  matrix[N,2] mu;
  
  // Sampling for Y
  mu[1,] = (mu_0' + X_centered[1,] * B_mat);
  // Y[1,] ~ multi_normal(mu[1,], Sigma);
  for (t in 2:N) {
    mu[t,] = (mu_0' + Y[t-1,] * auto_cor_mat + X_centered[t,] * B_mat);
    // Y[t,] ~ multi_normal(mu[t,], Sigma);
  }
  Y[,1] ~ normal(mu[,1] + sigma_w * rho_w * (Y[,2] - mu[,2]), sqrt(sigma_w^2 * (1 - rho_w^2)));
  Y[,2] ~ normal(mu[,2], 1);
  
  // Latent Length sampling
  matrix[N_obs, 2] mu_obs = mu[ind_obs,];
  vector[N_obs] lengths_llpd;
  for (t in 1:N_obs) {
      real A = quad_form(inverse(Sigma), U_obs[t,]');
      real B = U_obs[t,] * inverse(Sigma) * mu_obs[t,]';
      lengths_llpd[t] = log(latent_lengths[t]) - 1.0 / 2 * A * (latent_lengths[t] - B / A)^2;
  }
  target += sum(lengths_llpd);
}

generated quantities {
  real<lower=-pi(), upper=pi()> theta = atan2(mu_0[2], mu_0[1]);
  
  // Imputations for missing observations
  matrix[N_mis, 2] U_mis;
  U_mis = diag_pre_multiply(diagonal(Y_mis * Y_mis'), Y_mis);
  vector<lower=-pi(), upper=pi()>[N_mis] theta_imps;
  for (t in 1:N_mis) {
    theta_imps[t] = atan2(U_mis[t,2], U_mis[t,1]);
  }
  
  // Draws for (1,...,T) (Not one-shot prediction)
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
   vector<lower=-pi(), upper=pi()>[N] theta_ppd;
   theta_ppd[1] = atan2(U_ppd[1,2], U_ppd[1,1]);
   
  for (t in 2:(N+N_ppd)) {
    mu_ppd[t,] = (mu_0 + auto_cor_mat * Y[t-1,]')';
    // Prediction/backcast for time t given observed Y_1,...,Y_{t-1}
    Y_ppd[t, ] = multi_normal_rng(mu_ppd[t,] + eta[t,], Sigma)';
    U_ppd[t, ] = Y_ppd[t,] / sqrt(dot_self(Y_ppd[t,]));
    theta_ppd[t] = atan2(U_ppd[t, 2], U_ppd[t, 1]);
  }
}
