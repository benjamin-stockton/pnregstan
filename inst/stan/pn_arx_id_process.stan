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
}

transformed data {
  vector[K] X_bar = 1.0 / N * (rep_vector(1.0, N)' * X)';
  matrix[N, K] X_centered = X - rep_matrix(X_bar', N);
  
  // matrix[N, 2] U_lag = append_row(rep_vector(0,2), U[1:N,]);
  // X_centered = append_col(X_centered, U_lag);
  
  vector[K] X_ppd_bar = 1.0 / N_ppd * (rep_vector(1.0, N_ppd)' * X_ppd)';
  matrix[N_ppd, K] X_ppd_centered = X_ppd - rep_matrix(X_ppd_bar', N_ppd);
}

parameters {
  row_vector[2] mu_0;
  matrix[K,2] B_mat;
  matrix[2,2] auto_cor_mat;
  vector<lower=0>[N] latent_lengths;
}

transformed parameters {
  matrix[N, 2] Y;
  
  Y = diag_pre_multiply(latent_lengths, U);
}

model {
  matrix[N, 2] mu = rep_matrix(mu_0, N) + X_centered * B_mat;
  latent_lengths ~ normal(0, 10);
  
<<<<<<< Updated upstream
  mu_0 ~ normal(0, 1000);
  to_vector(B_mat) ~ normal(0, 1000);
  to_vector(auto_cor_mat) ~ normal(0, 1000);
=======
  mu_0 ~ normal(0, 10);
  to_vector(B_mat) ~ normal(0, 100);
  to_vector(auto_cor_mat) ~ normal(0, 1);
>>>>>>> Stashed changes
  
  // mu[1,] += mu_0;
  real A = dot_self(U[1,]);
  real B = U[1,] * mu[1,]';
<<<<<<< Updated upstream
  vector[N] lengths_llpd;
  lengths_llpd[1] = log(latent_lengths[1]) - 1.0 / 2 * A * (latent_lengths[1] - B / A)^2;
=======
  target += log(latent_lengths[1]) - 1.0 / 2 * A * (latent_lengths[1] - B / A)^2;
  
  // Y[1,] ~ multi_normal(mu[1,], identity_matrix(2));
  // Y[1,1] ~ normal(mu[1,1], 1);
  // Y[1,2] ~ normal(mu[1,2], 1);
>>>>>>> Stashed changes
  
  // Latent Length sampling
  for (t in 2:N) {
      mu[t,] += Y[t-1,] * auto_cor_mat;
      A = dot_self(U[t,]);
      B = U[t,] * mu[t,]';
      lengths_llpd[t] = log(latent_lengths[t]) - 1.0 / 2 * A * (latent_lengths[t] - B / A)^2;
      // Y[t,] ~ multi_normal_prec(mu[t,],  identity_matrix(2));
<<<<<<< Updated upstream
      // Y[t,1] ~ normal(mu[1], 1);
      // Y[t,2] ~ normal(mu[2], 1);
  }
  target += sum(lengths_llpd);
  
  // Y[1,] ~ multi_normal(mu[1,], identity_matrix(2));
  Y[,1] ~ normal(mu[,1], 1);
  Y[,2] ~ normal(mu[,2], 1);
=======
      
  }
    Y[,1] ~ normal(mu[,1], 1);
    Y[,2] ~ normal(mu[,2], 1);
>>>>>>> Stashed changes
}

generated quantities {// Post-processed Parameters
  
  real<lower=-pi(), upper=pi()> theta = atan2(mu_0[2], mu_0[1]);
  
  // Draws for (1,...,T,T+1,...,T+T_ppd) (Not one-shot prediction)
  matrix[N+N_ppd, K] X_mat = append_row(X_centered, X_ppd_centered);
  matrix[N+N_ppd, 2] mu_ppd;
  mu_ppd[1,] = mu_0;
  
  matrix[N+N_ppd,2] eta = X_mat * B_mat;
  
  matrix[N+N_ppd, 2] Y_ppd;
  // Y_ppd[1,] = multi_normal_rng(mu_ppd[1,] + eta[1,],  identity_matrix(2))';
  Y_ppd[1,1] = normal_rng(mu_ppd[1,1] + eta[1,1],  1);
  Y_ppd[1,2] = normal_rng(mu_ppd[1,2] + eta[1,2],  1);
  
  // ppd draws projected to S^1
  matrix[N+N_ppd, 2] U_ppd;
  U_ppd[1,] = Y_ppd[1,] / sqrt(dot_self(Y_ppd[1,]));
  
   // ppd draws converted to angles
   vector<lower=-pi(), upper=pi()>[N+N_ppd] theta_ppd;
   theta_ppd[1] = atan2(U_ppd[1,2], U_ppd[1,1]);
   
  for (t in 2:(N+1)) {
    mu_ppd[t,] = (mu_0 + Y[t-1,] * auto_cor_mat + eta[t,]);
    // Prediction/backcast for time t given observed Y_1,...,Y_{N+1}
    // Y_ppd[t, ] = multi_normal_rng(mu_ppd[t,] + eta[t,],  identity_matrix(2))';
    Y_ppd[t,1] = normal_rng(mu_ppd[t,1], 1);
    Y_ppd[t,2] = normal_rng(mu_ppd[t,2], 1);
    U_ppd[t, ] = Y_ppd[t,] / sqrt(dot_self(Y_ppd[t,]));
    theta_ppd[t] = atan2(U_ppd[t, 2], U_ppd[t, 1]);
  }
  
  for (t in 2:N_ppd) {
    mu_ppd[N+t,] = (mu_0 + Y_ppd[N+t-1,] * auto_cor_mat + eta[N+t,]);
    // Forecast for time t given observed Y_N+2,...,Y_{N+N_ppd}
    // Y_ppd[N+t, ] = multi_normal_rng(mu_ppd[N+t,] + eta[N+t,],  identity_matrix(2))';
    Y_ppd[N+t,1] = normal_rng(mu_ppd[N+t,1], 1);
    Y_ppd[N+t,2] = normal_rng(mu_ppd[N+t,2], 1);
    U_ppd[N+t, ] = Y_ppd[N+t,] / sqrt(dot_self(Y_ppd[N+t,]));
    theta_ppd[N+t] = atan2(U_ppd[N+t, 2], U_ppd[N+t, 1]);
  }
  
}
