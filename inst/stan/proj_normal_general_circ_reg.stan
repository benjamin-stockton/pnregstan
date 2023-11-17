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
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0>[N] latent_lengths;
  matrix[P, 2] B;
  vector[P] alpha;
  real<lower=0> sigma2;
  
  // real mu_B;
  // real<lower=0> tau_B;
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
    
    vector[N] gamma = X * alpha;
    
    // Priors
    matrix[N,2] mu = X * B;

    // mu_B ~ normal(0, 10);

    to_vector(B) ~ normal(0, 10000);
    alpha ~ normal(0, 10000);
    sigma2 ~ inv_gamma(0.01, 0.01);
    
    matrix[2, 2] Sigma_i = I_2;
    
    // Latent Length sampling
    for (i in 1:N) {
        Sigma_i[1, 1] = sigma2 + gamma[i]^2;
        Sigma_i[1, 2] = gamma[i];
        Sigma_i[2, 1] = gamma[i];
        Sigma_i[2, 2] = 1.0;
        matrix[2,2] Lambda_i = inverse_spd(Sigma_i);
        real A_i = quad_form(Lambda_i, U[i,]');
        real B_i = U[i,] * Lambda_i *  mu[i,]';
        target+= log(latent_lengths[i]) - 1.0 / 2 * A_i * (latent_lengths[i] - B_i / A_i)^2;
    }
    
    // matrix[P,P] Lambda_F = 1.0/ 5 * I_p + X' * X;
    
    // vector[P] mu_F_1 = inverse(Lambda_F) * X' * Y[,1];
    // vector[P] mu_F_2 = inverse(Lambda_F) * X' * Y[,2];
    
    // // Sampling for Y
    for (i in 1:N) {
      Sigma_i[1, 1] = sigma2 + gamma[i]^2;
      Sigma_i[1, 2] = gamma[i];
      Sigma_i[2, 1] = gamma[i];
      Sigma_i[2, 2] = 1;
      Y[i,] ~ multi_normal(mu[i,], Sigma_i);
    }
    // B[,1] ~ multi_normal_prec(mu_F_1, Lambda_F);
    // B[,2] ~ multi_normal_prec(mu_F_2, Lambda_F);
}

generated quantities {
    matrix[N, 2] mu = X * B;
    vector[N] gamma = X * alpha;
    vector[N] log_lik;
    matrix[N_tilde, 2] mu_tilde = X_tilde * B;
    vector[N_tilde] gamma_tilde = X_tilde * alpha;
    matrix[2, 2] Sigma_i = I_2;
    vector[N] A_i;
    vector[N] B_i;
    vector[N] C_i;
    vector[N] D_i;
    
    for (n in 1:N) {
        Sigma_i[1, 1] = sigma2 + gamma[n]^2;
        Sigma_i[1, 2] = gamma[n];
        Sigma_i[2, 1] = gamma[n];
        Sigma_i[2, 2] = 1;
        matrix[2,2] Lambda_i = inverse_spd(Sigma_i);
        
        A_i[n] = quad_form(Lambda_i, U[n,]');
        B_i[n] = U[n,] * Lambda_i *  mu[n,]';
        C_i[n] = -1.0/2 * quad_form(Lambda_i, mu[n,]');
        D_i[n] = B_i[n] / sqrt(A_i[n]);
        log_lik[n] = -log(2 * pi() * A_i[n]) - 1.0/2 * log(sigma2) + C_i[n] + log(1 + D_i[n] * exp(std_normal_lcdf(D_i[n]) - std_normal_lpdf(D_i[n])));
    }

    // for (n in 1:N_tilde)
    //     mu_tilde[n,] = X_tilde[n,] * B;

    // ppd draws from MVN
    matrix[N_tilde, 2] Y_tilde;
    // ppd draws projected to S^1
    matrix[N_tilde, 2] U_tilde;
     // ppd draws converted to angles
    vector<lower=-pi(), upper=pi()>[N_tilde] theta_tilde;

    for (n in 1:N_tilde) {
        Sigma_i[1, 1] = sigma2 + gamma_tilde[n]^2;
        Sigma_i[1, 2] = gamma_tilde[n];
        Sigma_i[2, 1] = gamma_tilde[n];
        Sigma_i[2, 2] = 1;
        Y_tilde[n,] = multi_normal_rng(mu_tilde[n,], Sigma_i)';
        U_tilde[n,] = Y_tilde[n,] / sqrt(dot_self(Y_tilde[n,]));
        theta_tilde[n] = atan2(U_tilde[n,2], U_tilde[n,1]);
    }
}


