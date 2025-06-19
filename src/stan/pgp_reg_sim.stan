data {
  int<lower=1> N;
  int<lower=0> P;
  array[N] vector[2] loc;
  matrix[N, P] X;
  real<lower=0> sigma_w;
  real rho_w;
  real rho;
  vector[2] mu0;
  matrix[P, 2] B;
}
transformed data {
  matrix[N, N] K = gp_exponential_cov(loc, 1.0, rho);
  matrix[N, N] L_K = cholesky_decompose(K);
  matrix[2,2] L_C = identity_matrix(2);
  L_C[1,1] = sigma_w;
  L_C[2,1] = rho_w;
  L_C[2,2] = sqrt(1 - rho_w^2);
  matrix[N,2] mu_mat = rep_matrix(mu0', N) + X * B;
  // for (n in 1:N) {
  //   K[n, n] = K[n, n];
  // }
}
parameters {
  matrix[N, 2] eta;
}
model {
  to_vector(eta) ~ std_normal();
}
generated quantities {
   matrix[N,2] y;
  y = mu_mat + L_K * eta * L_C';
  
  vector[N] theta;
  for (i in 1:N) {
    theta[i] = atan2(y[i,2], y[i,1]);
  }
}
