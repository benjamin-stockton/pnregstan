data {
  // Constants defining data dimensions
  int<lower=0> N_obs;       // number of data items
  int<lower=0> N_inc;       // number of incomplete cases
  int<lower=0> N_tilde;     // number of ppd draws
  int<lower=0> K;           // dimensionality of the directional covariates
  int<lower=0> P;           // number of predictors (not including directional ones)
  int<lower=1, upper = N_obs + N_inc> ii_obs[N_obs];
  int<lower=1, upper = N_obs + N_inc> ii_inc[N_inc];
  
  // Data matrices and vectors
  matrix[N_obs, P] X_obs;   // predictor matrix from the complete cases
  matrix[N_inc, P] X_inc;   // predictor matrix from the incomplete cases
                            //      not including the direcitonal covariates
  matrix[N_tilde, P] X_tilde; // predictor matrix for the ppd draws
  matrix[N_obs, K] U_obs;       // Observed directional covariates in cartesian coord
  matrix[N_tilde, K] U_tilde; // directional covariates for the ppd draws
  vector[N_obs] y_obs;      // outcome vector from the complete cases
  vector[N_inc] y_inc;      // outcomes from the incomplete cases
  
  // hyperparameters
  real<lower=0> tau_alpha;
  real<lower=0> tau_sigma;
  real<lower=0> tau_B;
  vector[P+K] mu_beta;
  real<lower=0> tau_beta;
}

transformed data {
  int<lower=0> N = N_obs + N_inc;
  matrix[N, P] X;
  vector[N] y;
  
  X[ii_obs,] = X_obs;
  X[ii_inc,] = X_inc;
  y[ii_obs] = y_obs;
  y[ii_inc] = y_inc;
  
  vector[K] ones;
  for (k in 1:K) {
      ones[k] = 1;
  }
  
  matrix[K, K] I_k = diag_matrix(ones);
}

parameters {
  real alpha;             // intercept
  vector[P+K] beta;       // coefficients for predictors
  real<lower=0> sigma;    // error scale
  // real<lower=0> sigma_alpha; // error for alpha
  // real<lower=0> sigma_beta; // error for beta
  
  vector<lower=0>[N_obs] latent_lengths;
  matrix[P+1, K] B;
  
  matrix[N_inc, K] V_inc;  // Imputed real vectors for missing angles; 
                           // Gets transformed to unit vectors in the next block
}

transformed parameters {
  matrix[N_obs, K] V_obs;
  matrix[N_inc, K] U_inc;
  
  V_obs = diag_pre_multiply(latent_lengths, U_obs);
  
  for (i in 1:N_inc) {
    U_inc[i,] = V_inc[i,] / sqrt(dot_self(V_inc[i,]));
  }
  
  matrix[N, K] U;
  
  U[ii_obs,] = U_obs;
  U[ii_inc,] = U_inc;
}

model {
  // Priors
  sigma ~ lognormal(0, tau_sigma);
  
  to_vector(B) ~ normal(0, tau_B);
  beta ~ normal(mu_beta, tau_beta);
  alpha ~ normal(0, tau_alpha);
  
  // Temporary mean transformations
  
  matrix[N_obs, K] eta_obs;
  matrix[N_inc, K] eta_inc;
  
  for (i in 1:N_obs) {
    eta_obs[i] = append_col(X_obs, y_obs)[i] * B;
  }
  for (i in 1:N_inc) {
    eta_inc[i] = append_col(X_inc, y_inc)[i] * B;
  }
  
  // Impute missing angles
  // Latent length sampling
  for (i in 1:N_obs) {
    real a = dot_self(U_obs[i,]);
    real b = dot_product(U_obs[i,], eta_obs[i,]);
    target += (K-1) * log(latent_lengths[i]) - 1.0/2 * a * (latent_lengths[i] - b / a)^2;
  }
  
  // Sampling for V_obs
  for (i in 1:N_obs) {
    V_obs[i,] ~ multi_normal(eta_obs[i], I_k);
  }
  // Sampling for V_inc
  for (i in 1:N_inc) {
    V_inc[i,] ~ multi_normal(eta_inc[i], I_k);
  }
  
  // Sampling for the response
  y ~ normal(alpha + append_col(X, U) * beta, sigma);  // likelihood
}

generated quantities {
  vector[N_tilde] y_tilde;
  vector[N_tilde] mu_tilde;
  
  mu_tilde = alpha + append_col(X_tilde, U_tilde) * beta;
  for (i in 1:N_tilde) {
    y_tilde[i] = normal_rng(mu_tilde[i], sigma);
  }
  
  vector[N_inc] theta_imp;
  
  for (i in 1:N_inc) {
    theta_imp[i] = atan2(U_inc[i, 2], U_inc[i, 1]);
  }
  
}

