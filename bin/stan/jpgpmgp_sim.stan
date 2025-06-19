data {
  int<lower=1> N;
  int<lower=0> P;
  array[N] vector[2] loc;
  real<lower=0> sigma_w;
  real<lower=-1,upper=1> rho_w;
  real<lower=0> rho;
  
  vector<lower=0>[P] sigma_x;
  vector<lower=0>[P] alpha_x;
  real<lower=0> rho_x;
  vector[2] mu_w;
  matrix[3, P] B_x;
}
transformed data {
  matrix[N, N] K = gp_exponential_cov(loc, 1.0, rho);
  matrix[N, N] L_K = cholesky_decompose(K);
  matrix[2,2] L_C = identity_matrix(2);
  L_C[1,1] = sigma_w;
  L_C[2,1] = rho_w;
  L_C[2,2] = sqrt(1 - rho_w^2);
  matrix[N,2] mu_mat = rep_matrix(mu_w', N);
  
  matrix[P, P] Sigma_x = diag_matrix(sigma_x^2);
  
  matrix[N, N] K_x = gp_exponential_cov(loc, 1.0, rho_x);
  matrix[N, N] L_K_x = cholesky_decompose(K);
  matrix[P,P] L_Sigma = cholesky_decompose(Sigma_x);
  // for (n in 1:N) {
  //   K[n, n] = K[n, n];
  // }
}
parameters {
  matrix[N, 2] eta;
  
  matrix[N, P] nu;
}
model {
  to_vector(eta) ~ std_normal();
  to_vector(nu) ~ std_normal();
}
generated quantities {
   matrix[N,2] y;
  y = mu_mat + L_K * eta * L_C';
  
  vector[N] theta;
  for (i in 1:N) {
    y[i,] = multi_normal_rng(y[i,]', tcrossprod(L_C))';
    theta[i] = atan2(y[i,2], y[i,1]);
  }
  
  matrix[N,P] mu_mat_x;
  for (i in 1:N) {
    mu_mat_x[i,] = B_x[1,] + y[i,] * B_x[2:3,];
  }
  matrix[N,P] x;
  x = mu_mat_x + L_K_x * nu * diag_pre_multiply(alpha_x, L_Sigma)';
}
