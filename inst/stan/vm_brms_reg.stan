// generated with brms 2.18.0
functions {
  /* compute the tan_half link
   * Args:
   *   x: a vector in (-pi, pi)
   * Returns:
   *   a vector in (-Inf, Inf)
   */
   vector tan_half(vector x) {
     return tan(x / 2);
   }
  /* compute the inverse of the tan_half link
   * Args:
   *   y: a vector in (-Inf, Inf)
   * Returns:
   *   a vector in (-pi, pi)
   */
   vector inv_tan_half(vector y) {
     return 2 * atan(y);
   }
  /* von Mises log-PDF of a single response
   * for kappa > 100 the normal approximation is used
   * for reasons of numerial stability
   * Args:
   *   y: the response vector between -pi and pi
   *   mu: location parameter vector
   *   kappa: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real von_mises_real_lpdf(real y, real mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(y | mu, kappa);
     } else {
       return normal_lpdf(y | mu, sqrt(1 / kappa));
     }
   }
  /* von Mises log-PDF of a response vector
   * for kappa > 100 the normal approximation is used
   * for reasons of numerial stability
   * Args:
   *   y: the response vector between -pi and pi
   *   mu: location parameter vector
   *   kappa: precision parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real von_mises_vector_lpdf(vector y, vector mu, real kappa) {
     if (kappa < 100) {
       return von_mises_lpdf(y | mu, kappa);
     } else {
       return normal_lpdf(y | mu, sqrt(1 / kappa));
     }
   }
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=0> N_ppd;
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  matrix[N_ppd, K] X_ppd;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  matrix[N_ppd, Kc] X_ppd_c;
  vector[Kc] means_X;  // column means of X before centering
  vector[Kc] means_X_ppd;
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 1:Kc) {
    means_X_ppd[i] = mean(X_ppd[, i]);
    X_ppd_c[, i] = X_ppd[, i] - means_X_ppd[i];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> kappa;  // precision parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b | 0, 5);
  lprior += normal_lpdf(Intercept | 0, 2);
  lprior += gamma_lpdf(kappa | 2, 0.01);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Xc * b;
    mu = Intercept + inv_tan_half(mu);
    target += von_mises_vector_lpdf(Y | mu, kappa);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  
  vector[N_ppd] mu_ppd = rep_vector(0.0, N_ppd);
  mu_ppd += X_ppd_c * b;
  mu_ppd = Intercept + inv_tan_half(mu_ppd);
  
  vector[N_ppd] theta_ppd;
  for (i in 1:N_ppd) {
     theta_ppd[i] = von_mises_rng(mu_ppd[i], kappa);
  }
}

