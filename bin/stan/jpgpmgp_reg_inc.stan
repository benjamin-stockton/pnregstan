functions {

  matrix build_cholesky_factor(vector y1,
                                 array[] vector loc1,
                                 real alpha,
                                 real rho,
                                 real sigma) {
        int N1 = rows(y1);
        
        matrix[N1, N1] K = gp_exponential_cov(loc1, alpha, rho);
        real sq_sigma = square(sigma);
        
        K = add_diag(K, sq_sigma);
        
        matrix[N1, N1] L_K = cholesky_decompose(K);
        return L_K;
    }


  vector gp_pred_rng(array[] vector loc2,
                     vector y1,
                     array[] vector loc1,
                     vector mu1,
                     vector mu2,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = rows(y1);
    int N2 = size(loc2);
    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] y_minus_mu;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      K = gp_exponential_cov(loc1, alpha, rho);
      for (n in 1:N1) {
        K[n, n] = K[n, n] + square(sigma);
      }
      L_K = cholesky_decompose(K);
      y_minus_mu = y1 - mu1;
      K_div_y1 = mdivide_left_tri_low(L_K, y_minus_mu);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = gp_exponential_cov(loc1, loc2, alpha, rho);
      f2_mu = mu2 + (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = gp_exponential_cov(loc2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }


  matrix kronecker_prod(matrix A, matrix B) {
      matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
      int m;
      int n;
      int p;
      int q;
      m = rows(A);
      n = cols(A);
      p = rows(B);
      q = cols(B);
      for (i in 1:m) {
        for (j in 1:n) {
          int row_start;
          int row_end;
          int col_start;
          int col_end;
          row_start = (i - 1) * p + 1;
          row_end = (i - 1) * p + p;
          col_start = (j - 1) * q + 1;
          col_end = (j - 1) * q + q;
          C[row_start:row_end, col_start:col_end] = A[i, j] * B;
        }
      }
      return C;
    }

vector latent_length_llpd(vector latent_lengths,
                          matrix W1,
                          matrix mu_mat,
                          matrix K_inv,
                          matrix Sigma_inv,
                          int N1) {
  vector[N1] lengths_llpd;
  real a;
  real b;
  vector[2] E_si;
  real w_tmp;
  for (i in 1:N1) {
    vector[2] U1i = W1[i,]';
    U1i = U1i / sqrt(U1i[1]^2 + U1i[2]^2);
    w_tmp = 1/K_inv[i,i];
    vector[2] mu_w = mu_mat[i,]';
    E_si = mu_w - w_tmp * (W1 - mu_mat)' * K_inv[,i] + w_tmp * (W1[i,] - mu_w')' * K_inv[i,i];

    a = quad_form(Sigma_inv, U1i);
    b = U1i' * Sigma_inv * E_si;
    lengths_llpd[i] = log(latent_lengths[i]) - K_inv[i,i] / 2.0 * (a * latent_lengths[i] - b)^2;
    // lengths_llpd[i] = log(latent_lengths[i]) - K_inv[i,i] / 2.0 * (latent_lengths[i]^2 * a - 2 * latent_lengths[i] * b);

  }
  return lengths_llpd;
}

}

data {
 int<lower=1> N_obs; //number of data points
 int<lower=0> N_inc;
 int<lower=1> N2;
 int<lower=0> P; // number of inline predictors
 vector[N_obs] theta_obs; // angular outcomes
 vector[N_inc] theta_inc;
 matrix[N_obs, P] Xobs;
 matrix[N_inc, P] Xinc;
 vector[N_obs] Yobs;
 array[N_obs] vector[2] loc_obs;
 array[N_inc] vector[2] loc_inc;
 array[N2] vector[2] loc2;
}

transformed data{
 int<lower=1> N;
 N = N_inc + N_obs + N2;
 int N1 = N_inc + N_obs;
 real delta = 1e-9;
 matrix[N1, 2] U1;
 for (i in 1:N_obs) {
   U1[i,1] = cos(theta_obs[i]);
   U1[i,2] = sin(theta_obs[i]);
 }
 for (i in 1:N_inc) {
   U1[N_obs + i,1] = cos(theta_inc[i]);
   U1[N_obs + i,2] = sin(theta_inc[i]);
 }
 
 array[N] vector[2] loc;
 for (n in 1:N_obs) {
   loc[n] = loc_obs[n];
 }
 for (n in 1:N_inc) {
   loc[N_obs+n] = loc_inc[n];
 }
 for (n in 1:N2) {
   loc[N1+n] = loc2[n];
 }
 // matrix[2,2] L_Sigma = identity_matrix(2);
 // vector[2] mu0 = rep_vector(0, 2);
 vector[2] alpha = rep_vector(1.0, 2);
 // matrix[P,2] B = rep_matrix(rep_vector(0, 2)', P);
 // real rho = 0.91;
 vector[N1] zeros = rep_vector(0, N1);
 
 matrix[N1, P] X1 = append_row(Xobs, Xinc);
}

parameters{
 vector<lower=0>[N1] latent_lengths;
 vector[2] mu0;
 
 real<lower=0> rho;
 // vector<lower=0>[2] alpha;
 // cholesky_factor_corr[2] L_Sigma;
 matrix[N1, 2] eta;
 
 matrix[P, 2] B;
 
 real<lower=0> sigma_w;
 real<lower=-1, upper=1> rho_w;
 
 vector[P] mu_x;
 vector<lower=0>[P] sigma_x;
 vector<lower=0>[P]  rho_x;
 vector<lower=0>[P] alpha_x;
 
 
 vector[P+3] beta_y;
 real<lower=0> sigma_y;
 real<lower=0> rho_y;
 real<lower=0> alpha_y;
 
 vector[N_inc] Yinc;
}

transformed parameters {
  matrix[N1, 2] W1;

  W1 = diag_pre_multiply(latent_lengths, U1);
 
 matrix[N1, 2] mu_mat_w;
  for (i in 1:N1) {
    mu_mat_w[i,] = mu0' + X1[i,] * B;
  }
  
  real sigma_w_sq = sigma_w^2;
  
  matrix[2,2] L_Sigma;
  L_Sigma[1,1] = sigma_w;
  L_Sigma[1,2] = 0;
  L_Sigma[2,1] = rho_w;
  L_Sigma[2,2] = sqrt(1.0 - rho_w^2);
  
  vector[N1] Y1 = append_row(Yobs, Yinc);
}

model{
    mu_x ~ normal(0, 5);
    sigma_x ~ inv_gamma(2,1);
    rho_x ~ inv_gamma(2,1);
    alpha_x ~ normal(0,5);

    for (p in 1:P) {
      matrix[N1, N1] L_K = build_cholesky_factor(X1[,p], loc[1:N1,], alpha_x[p], rho_x[p], sigma_x[p]);
      vector[N1] mu_p = rep_vector(mu_x[p], N1);
      
      X1[,p] ~ multi_normal_cholesky(mu_p, L_K);
    }
    
    
    mu0 ~ normal(0, 10);
    to_vector(B) ~ normal(0, 10);
    rho ~ inv_gamma(2,1);
    sigma_w_sq ~ inv_gamma(2,1);
    // rho_w ~ std_normal();
    to_vector(eta) ~ std_normal();
    latent_lengths ~ normal(0, 10);
    
    matrix[N1, 2] f;
    matrix[N1, N1] K;
    matrix[N1, N1] K_inv;
    matrix[2,2] Sigma_inv;
    Sigma_inv = chol2inv(L_Sigma);
    {
      K = gp_exponential_cov(loc[1:N1,], 1.0, rho);
      matrix[N1, N1] L_K;

      // diagonal elements
      for (n in 1:N1) {
        K[n, n] = K[n, n] + delta;
      }

      L_K = cholesky_decompose(K);
      K_inv = chol2inv(L_K);
      f = L_K * eta * L_Sigma';
    }
    
    // Latent Length sampling
    vector[N1] lengths_llpd;
    {
      real a;
      real b;
      vector[2] E_si;
      for (i in 1:N1) {
        real w_tmp;
        w_tmp = 1/K_inv[i,i];
        E_si = mu_mat_w[i,]' - w_tmp * (W1 - mu_mat_w[1:N1,])' * K_inv[1:N1,i] + w_tmp * (W1[i,] - mu_mat_w[i,])' * K_inv[i,i];
      
        a = quad_form(Sigma_inv, U1[i,]');
        b = U1[i,] * Sigma_inv * E_si;
        lengths_llpd[i] = log(latent_lengths[i]) - K_inv[i,i] / 2.0 * a * (latent_lengths[i] - b / a)^2;
      }
    }
    target += sum(lengths_llpd);
    
    to_vector(W1) ~ normal(to_vector(mu_mat_w + f), 1.0);
    
    
    beta_y ~ normal(0, 5);
    sigma_y ~ inv_gamma(2,1);
    rho_y ~ inv_gamma(2,1);
    alpha_y ~ normal(0,5);
    
    matrix[N1, N1] L_K_y = build_cholesky_factor(Y1, loc[1:N1,], alpha_y, rho_y, sigma_y);
    
    vector[N1] mu_y = beta_y[1] + X1 * beta_y[2:(P+1)] + W1 * beta_y[(P+2):(P+3)];

    Y1 ~ multi_normal_cholesky(mu_y, L_K_y);
}

generated quantities{
    real spat_decay = 1.0 / rho;
    
    // In-Sample PPD Draws
    matrix[N1, P] x_pred1;
    {
      
    for (p in 1:P) {
      vector[N1] mu_pred_0 = rep_vector(mu_x[p], N1);
      vector[N1] f1_x;
      f1_x = gp_pred_rng(loc[1:N1,],
                     X1[,p],
                     loc[1:N1,],
                     mu_pred_0,
                     mu_pred_0,
                     alpha_x[p],
                     rho_x[p],
                     sigma_x[p],
                     delta);

      for (n in 1:N1) {
        x_pred1[n,p] = normal_rng(f1_x[n], sigma_x[p]);
      }

    }
    }

    matrix[N1, 2] w_pred1;
    vector[N1] theta_pred1;
    {
      matrix[N1, 2] f_1;
      matrix[N1, N1] K = gp_exponential_cov(loc[1:N1,], 1.0, rho);
      matrix[N1, N1] L_K;

      // diagonal elements
      for (n in 1:N1) {
        K[n, n] = K[n, n] + delta;
      }

      L_K = cholesky_decompose(K);
      f_1 = L_K * eta * L_Sigma';

      for(i in 1:N1) {
          // w_pred1[i,] = multi_normal_rng(mu_mat_w[i,]' + f_1[i,]', identity_matrix(2))';
          w_pred1[i,] = multi_normal_rng(mu_mat_w[i,]' + f_1[i,]', tcrossprod(L_Sigma))';
          theta_pred1[i] = atan2(w_pred1[i,2], w_pred1[i,1]);
      }
    }
    
    vector[N1] y_pred1;
    {
      vector[N1] mu_pred_y = beta_y[1] + X1 * beta_y[2:(P+1)] + W1 * beta_y[(P+2):(P+3)];
      vector[N1] f1_y;
      f1_y = gp_pred_rng(loc[1:N1,],
                     Y1,
                     loc[1:N1,],
                     mu_pred_y,
                     mu_pred_y,
                     alpha_y,
                     rho_y,
                     sigma_y,
                     delta);

      for (n in 1:N1) {
        y_pred1[n] = normal_rng(f1_y[n], sigma_y);
      }
    }

    // Regular Grid/Out-of-Sample PPD Draws
    matrix[N2,P] x_pred2;
    {
    for (p in 1:P) {
      vector[N2] f2_x;
      vector[N1] mu_pred_0 = rep_vector(mu_x[p], N1);
      vector[N2] mu_pred_02 = rep_vector(mu_x[p], N2);
      f2_x = gp_pred_rng(loc2,
                     X1[,p],
                     loc[1:N1,],
                     mu_pred_0,
                     mu_pred_02,
                     alpha_x[p],
                     rho_x[p],
                     sigma_x[p],
                     delta);
      for (n in 1:N2) {
        x_pred2[n,p] = normal_rng(f2_x[n], sigma_x[p]);
      }
    }
    }

    matrix[N2, 2] w_pred2;
    vector[N2] theta_pred2;
    {
      vector[2*N2] vec_Y2;
      vector[2*N2] xi;
      matrix[N1, N1] K = gp_exponential_cov(loc[1:N1,], 1.0, rho);
      matrix[N1, N1] L_K;
      // diagonal elements
      for (n in 1:N1) {
        K[n, n] = K[n, n] + delta;
      }
      L_K = cholesky_decompose(K);

      matrix[N2, N2] K_pred = gp_exponential_cov(loc2, 1.0, rho);
      for (n in 1:N2) {
        K_pred[n,n] += delta;

        // draws from standard normal
        xi[n] = std_normal_rng();
        xi[N2+n] = std_normal_rng();
      }
      matrix[N2, N2] L_Kpred = cholesky_decompose(K_pred);

      matrix[N1, N2] K_loc1_loc2 = gp_exponential_cov(loc[1:N1,], loc2, 1.0, rho);

      vector[2*N1] v1_minus_mu = to_vector(W1' - mu_mat_w');
      matrix[N1, N2] U = mdivide_left_tri_low(L_K, K_loc1_loc2);

      matrix[N2, 2] mu2 = rep_matrix(mu0', N2) + x_pred2 * B;
      vector[2*N2] cond_mean;
      matrix[2*N2, 2*N2] cond_cov;

      cond_mean = to_vector(mu2') + kronecker_prod(mdivide_right_tri_low(U', L_K), identity_matrix(2)) * v1_minus_mu;

      matrix[2,2] Sigma = L_Sigma*L_Sigma';
      cond_cov = kronecker_prod(K_pred - crossprod(U), Sigma);

      vec_Y2 = cond_mean + cholesky_decompose(cond_cov) * xi;
      w_pred2 = to_matrix(vec_Y2, 2, N2)';
    }


    for(i in 1:N2) {
        theta_pred2[i] = atan2(w_pred2[i,2], w_pred2[i,1]);
    }
    
    vector[N2] y_pred2;
    {
    {
      vector[N2] f2_y;
      vector[N1] mu_pred_y = beta_y[1] + X1 * beta_y[2:(P+1)] + W1 * beta_y[(P+2):(P+3)];
      vector[N2] mu_pred_y2 = beta_y[1] + x_pred2 * beta_y[2:(P+1)] + w_pred2 * beta_y[(P+2):(P+3)];
      f2_y = gp_pred_rng(loc2,
                     Y1,
                     loc[1:N1,],
                     mu_pred_y,
                     mu_pred_y2,
                     alpha_y,
                     rho_y,
                     sigma_y,
                     delta);
      for (n in 1:N2) {
        y_pred2[n] = normal_rng(f2_y[n], sigma_y);
      }
    }
    }
}


