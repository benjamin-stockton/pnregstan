// The input data is a array 'U' of size 'N x 2' and covariate matrix
// X of size 'N X P'. 
data {
  int<lower=0> N;
  int<lower=0> N_ppd;
  int<lower=0> P;
  matrix[N,2] U;
  matrix[N,P] X;
  matrix[N_ppd, P] X_ppd;
}

transformed data {
  vector[P] X_bar = 1.0 / N * (rep_vector(1.0, N)' * X)';
  matrix[N, P] X_centered = X - rep_matrix(X_bar', N);
  X_centered[,1] = rep_vector(1, N);
  
  vector[P] X_ppd_bar = 1.0 / N_ppd * (rep_vector(1.0, N_ppd)' * X_ppd)';
  matrix[N_ppd, P] X_ppd_centered = X_ppd - rep_matrix(X_ppd_bar', N_ppd);
  X_ppd_centered[,1] = rep_vector(1, N_ppd);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0>[N] latent_lengths;
  matrix[P, 2] B;
  
  // real mu_B;
  // real<lower=0> tau_B;
}

transformed parameters {
    matrix[N,2] Y;
    
    Y = diag_pre_multiply(latent_lengths, U);
}

// The model to be estimated. We model the output
// 'Y' to be multivariate normal-ly distributed with mean 'mu'
// and identity covariance.
model {
    
    // Priors
    matrix[N,2] mu = X_centered * B;
    
    latent_lengths ~ normal(0, 10);

    // mu_B ~ normal(0, 10);
    // tau_B ~ inv_gamma(0.001, 0.001);

    to_vector(B) ~ normal(0, 10);
    
    // Latent Length sampling
    // vector[N] lengths_llpd;
    // vector[N] a = rows_dot_self(U);
    // vector[N] b = rows_dot_product(U, mu);
    for (i in 1:N) {
        real a = dot_self(U[i,]);
        real b = dot_product(U[i,], mu[i,]);
        // lengths_llpd[i] = log(latent_lengths[i]) - 1.0 / 2 * a * (latent_lengths[i] - b / a)^2;
        target += log(latent_lengths[i]) - 1.0 / 2 * (latent_lengths[i] * a - b)^2;
        Y[i,] ~ multi_normal(mu[i,], identity_matrix(2));
    }
    // target += sum(log(latent_lengths)) - 1.0 / 2 * quad_form(identity_matrix(N), (diag_matrix(a) * latent_lengths - b));
    // target += sum(lengths_llpd);
    
    // matrix[P,P] Lambda_F = 1.0/ 100 * identity_matrix(P) + quad_form(identity_matrix(N), X);
    // 
    // vector[P] mu_F_1 = inverse(Lambda_F) * X' * Y[,1];
    // vector[P] mu_F_2 = inverse(Lambda_F) * X' * Y[,2];
    // B[,1] ~ multi_normal_prec(mu_F_1, Lambda_F);
    // B[,2] ~ multi_normal_prec(mu_F_2, Lambda_F);
    
    // // Sampling for Y
    // Y ~ multi_normal(mu, I_2);
    // for (i in 1:N) {
    //   Y[i,] ~ multi_normal(mu[i,], identity_matrix(2));
    // }
}

generated quantities {
    matrix[N, 2] mu = X_centered * B;
    matrix[N_ppd, 2] mu_ppd = X_ppd_centered * B;
    // vector[N] log_lik;
    // 
    // for (i in 1:N) {
    //     real A_i = dot_self(U[i,]');
    //     real B_i = dot_product(U[i,]',  mu[i,]);
    //     real C = -1.0/2 * dot_self(mu[i,]');
    //     real D_i = B_i / sqrt(A_i);
    //     log_lik[i] = -log(2 * pi() * A_i) + C + log(1 + D_i * std_normal_cdf(D_i) / exp(std_normal_lpdf(D_i)));
    // }

    // ppd draws from MVN
    matrix[N_ppd, 2] Y_ppd;
    
    // Y_ppd[,1] = std_normal_rng() + mu_ppd[,1];
    // Y_ppd[,2] = std_normal_rng() + mu_ppd[,2];
    // ppd draws projected to S^1
    // vector[N_ppd] lens = rows_dot_self(Y_ppd);
    matrix[N_ppd, 2] U_ppd;
    // U_ppd = diag_matrix(1.0 / lens) * Y_ppd;
     // ppd draws converted to angles
    vector<lower=-pi(), upper=pi()>[N_ppd] theta_ppd;
    // theta_ppd = atan2(U_ppd[,2], U_ppd[,1]);
    
    for (n in 1:N_ppd) {
        Y_ppd[n,] = multi_normal_rng(mu_ppd[n,], identity_matrix(2))';
        U_ppd[n,] = Y_ppd[n,] / sqrt(dot_self(Y_ppd[n,]));
        theta_ppd[n] = atan2(U_ppd[n,2], U_ppd[n,1]);
    }
}
