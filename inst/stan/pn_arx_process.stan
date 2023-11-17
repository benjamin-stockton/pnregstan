//
// This Stan program fits the projected normal AR(1) model

data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N, 2] U;
  matrix[N, K] X;
  real<lower=0> sigma_0;
  
  int<lower=0> N_tilde;
  matrix[N_tilde, K] X_tilde;
}

parameters {
  vector[2] mu_0;
  matrix[K,2] B_mat;
  matrix[2,2] auto_cor_mat;
  vector<lower=0>[N] latent_lengths;
  
  real<lower=0> sigma_w;
  real<lower=-1, upper=1> rho_w;
}

transformed parameters {
  matrix[N, 2] Y;
  
  Y = diag_pre_multiply(latent_lengths, U);
  
  matrix[2,2] Sigma;
  Sigma[1,1] = sigma_w^2;
  Sigma[1,2] = sigma_w * rho_w;
  Sigma[2,1] = sigma_w * rho_w;
  Sigma[2,2] = 1;
}

model {
  matrix[N,2] mu;
  
  mu[1,] = (mu_0' + X[1,] * B_mat);
  for (t in 2:N) {
    mu[t,] = (mu_0' + Y[t-1,] * auto_cor_mat + X[t,] * B_mat);
  }
  
  // Latent Length sampling
  for (t in 1:N) {
      real A = quad_form(inverse_spd(Sigma), U[t,]');
      real B = U[t,] * inverse_spd(Sigma) * mu[t,]';
      target += log(latent_lengths[t]) - 1.0 / 2 * A * (latent_lengths[t] - B / A)^2;
  }
  
  // Sampling for Y
  for (t in 1:N) {
      Y[t,] ~ multi_normal(mu[t,], Sigma);
  }
}

generated quantities {
  real<lower=-pi(), upper=pi()> theta = atan2(mu_0[2], mu_0[1]);
  
  // Draws for (1,...,T) (Not one-shot prediction)
  matrix[N, 2] mu_ppd;
  mu_ppd[1,] = mu_0';
  
  matrix[N,2] eta = X * B_mat;
  
  matrix[N, 2] Y_ppd;
  Y_ppd[1,] = multi_normal_rng(mu_ppd[1,] + eta[1,], Sigma)';
  
  // ppd draws projected to S^1
  matrix[N, 2] U_ppd;
  U_ppd[1,] = Y_ppd[1,] / sqrt(dot_self(Y_ppd[1,]));
  
   // ppd draws converted to angles
   vector<lower=-pi(), upper=pi()>[N] theta_ppd;
   theta_ppd[1] = atan2(U_ppd[1,2], U_ppd[1,1]);
   
  for (t in 2:N) {
    mu_ppd[t,] = (mu_0 + auto_cor_mat * Y[t-1,]')';
    // Prediction/backcast for time t given observed Y_1,...,Y_{t-1}
    Y_ppd[t, ] = multi_normal_rng(mu_ppd[t,] + eta[t,], Sigma)';
    U_ppd[t, ] = Y_ppd[t,] / sqrt(dot_self(Y_ppd[t,]));
    theta_ppd[t] = atan2(U_ppd[t, 2], U_ppd[t, 1]);
  }
    
  matrix[N_tilde, 2] mu_tilde;
  mu_tilde[1,] = (mu_0' + Y[N,] * auto_cor_mat + X_tilde[1,] * B_mat);
    
  // ppd draws from MVN
  matrix[N_tilde, 2] Y_tilde;
  Y_tilde[1,] = multi_normal_rng(mu_tilde[1,], Sigma)';
  
  // ppd draws projected to S^1
  matrix[N_tilde, 2] U_tilde;
  U_tilde[1,] = Y_tilde[1,] / sqrt(dot_self(Y_tilde[1,]));
  
   // ppd draws converted to angles
   vector<lower=-pi(), upper=pi()>[N_tilde] theta_tilde;
   theta_tilde[1] = atan2(U_tilde[1,2], U_tilde[1,1]);
   
  for (t in 2:N_tilde) {
    mu_tilde[t,] = (mu_0' + Y_tilde[t-1,] * auto_cor_mat + X_tilde[t,] * B_mat);
    Y_tilde[t, ] = multi_normal_rng(mu_tilde[t,], Sigma)';
    U_tilde[t, ] = Y_tilde[t,] / sqrt(dot_self(Y_tilde[t,]));
    theta_tilde[t] = atan2(U_tilde[t, 2], U_tilde[t, 1]);
  }
}
