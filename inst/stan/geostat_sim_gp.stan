data {
  int<lower=1> N;
  array[N] vector[2] loc;
  real<lower=0> sigma;
  real<lower=0> alpha;
  real rho;
}
transformed data {
  matrix[N, N] K = gp_exponential_cov(loc, alpha, rho);
  vector[N] mu = rep_vector(0, N);
  for (n in 1:N) {
    K[n, n] = K[n, n] + sigma^2;
  }
  matrix[N, N] L_K = cholesky_decompose(K);
}
parameters {
  vector[N] eta;
}
model {
  eta ~ std_normal();
}
generated quantities {
  vector[N] y;
  y = mu + L_K * eta;
}
