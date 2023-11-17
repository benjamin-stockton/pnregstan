// The input data is a array 'U' of size 'N x 2' and covariate matrix
// X of size 'N X P'. 
data {
  int<lower=0> N;
  int<lower=0> N_tilde;
  int<lower=0> P;
  matrix[N,2] U;
  matrix[N,P] X;
  matrix[N_tilde, P] X_tilde;
}

transformed data {
    vector[2] ones;
    for (k in 1:2) {
        ones[k] = 1;
    }
    matrix[2, 2] I_2 = diag_matrix(ones);
    
    vector[P] ones2;
    for (k in 1:P) {
        ones2[k] = 1;
    }
    matrix[P, P] I_p = diag_matrix(ones2);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0>[N] latent_lengths;
  matrix[P, 2] B;
  
  real mu_B;
  real<lower=0> tau_B;
}

transformed parameters {
    matrix[N,2] Y;
    
    // Y = diag_pre_multiply(latent_lengths, U);
    Y = diag_matrix(latent_lengths) * U;
}

// The model to be estimated. We model the output
// 'Y' to be multivariate normal-ly distributed with mean 'mu'
// and identity covariance.
model {
    
    // Priors
    matrix[N,2] mu = X * B;

    mu_B ~ normal(0, 10);
    tau_B ~ inv_gamma(0.001, 0.001);

    to_vector(B) ~ normal(mu_B, tau_B);
    
    // Latent Length sampling
    for (i in 1:N) {
        real a = dot_self(U[i,]);
        real b = dot_product(U[i,], mu[i,]);
        target+= (2-1) * log(latent_lengths[i]) - 1.0 / 2 * a * (latent_lengths[i] - b / a)^2;
    }
    
    // matrix[P,P] Lambda_F = 1.0/ 5 * I_p + X' * X;
    
    // vector[P] mu_F_1 = inverse(Lambda_F) * X' * Y[,1];
    // vector[P] mu_F_2 = inverse(Lambda_F) * X' * Y[,2];
    
    // // Sampling for Y
    // Y ~ multi_normal(mu, I_2);
    for (i in 1:N) {
      Y[i,] ~ multi_normal(mu[i,], I_2);
    }
    // B[,1] ~ multi_normal_prec(mu_F_1, Lambda_F);
    // B[,2] ~ multi_normal_prec(mu_F_2, Lambda_F);
}

generated quantities {
    matrix[N, 2] mu = X * B;
    matrix[N_tilde, 2] mu_tilde = X_tilde * B;
    vector[N] log_lik;
    
    for (i in 1:N) {
        real A_i = dot_self(U[i,]');
        real B_i = dot_product(U[i,]',  mu[i,]);
        real C = -1.0/2 * dot_self(mu[i,]');
        real D_i = B_i / sqrt(A_i);
        log_lik[i] = -log(2 * pi() * A_i) + C + log(1 + D_i * std_normal_cdf(D_i) / exp(std_normal_lpdf(D_i)));
    }

    // ppd draws from MVN
    matrix[N_tilde, 2] Y_tilde;
    // ppd draws projected to S^1
    matrix[N_tilde, 2] U_tilde;
     // ppd draws converted to angles
    vector<lower=-pi(), upper=pi()>[N_tilde] theta_tilde;
    
    for (n in 1:N_tilde) {
        Y_tilde[n,] = multi_normal_rng(mu_tilde[n,], I_2)';
        U_tilde[n,] = Y_tilde[n,] / sqrt(dot_self(Y_tilde[n,]));
        theta_tilde[n] = atan2(U_tilde[n,2], U_tilde[n,1]);
    }
}
