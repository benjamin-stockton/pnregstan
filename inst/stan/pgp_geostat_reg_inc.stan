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
 int<lower=0> P;
 vector[N_obs] theta1; // angular outcomes
 matrix[N_obs, P] Xobs;
 matrix[N_inc, P] Xinc;
 array[N_obs] vector[2] loc_obs;
 array[N_inc] vector[2] loc_inc;
 array[N2] vector[2] loc2;
}

transformed data{
 int<lower=1> N;
 N = N_inc + N_obs + N2;
 int N1 = N_inc + N_obs;
 real delta = 1e-9;
 matrix[N_obs, 2] Uobs;
 for (i in 1:N_obs) {
   Uobs[i,1] = cos(theta1[i]);
   Uobs[i,2] = sin(theta1[i]);
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

 // matrix[2*N1, 2*P] X_aug = kronecker_prod(identity_matrix(2), X1);
 
 matrix[N1, P] Q_ast;
  matrix[P, P] R_ast;
  matrix[P, P] R_ast_inverse;
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X1) * sqrt(N1 - 1);
  R_ast = qr_thin_R(X1) / sqrt(N1 - 1);
  R_ast_inverse = inverse(R_ast);
 
 // hyperparameters
 real<lower=0> lambda0 = 100;
 real<lower=0> a_s = 2;
 real<lower=0> b_s = 1;
}

parameters{
 vector<lower=0>[N_obs] latent_lengths;
 vector[2] mu0;
 
 real<lower=0> rho;
 // vector<lower=0>[2] alpha;
 // cholesky_factor_corr[2] L_Sigma;
 matrix[N1, 2] eta;
 vector[2*N2] xi;
 
 matrix[P, 2] coef_mat;
 
 real<lower=0> sigma_w;
 real<lower=-1, upper=1> rho_w;
 
 vector[P] mu_x;
 vector<lower=0>[P] sigma_x;
 vector<lower=0>[P]  rho_x;
 vector<lower=0>[P] alpha_x;
 
 matrix[N_inc, 2] Winc;
}

transformed parameters {
  matrix[N_obs, 2] Wobs;

  Wobs = diag_pre_multiply(latent_lengths, Uobs);
  matrix[N1, 2] W1 = append_row(Wobs, Winc);
 
  matrix[N1, 2] mu_mat;
  mu_mat = rep_matrix(mu0', N1) + Q_ast * coef_mat;
  // for (i in 1:N1) {
  //   mu_mat[i,] = mu0' + X1[i,] * coef_mat;
  // }
  
  real sigma_w_sq = sigma_w^2;
  
  matrix[2,2] L_Sigma;
  L_Sigma[1,1] = sigma_w;
  L_Sigma[1,2] = 0;
  L_Sigma[2,1] = rho_w;
  L_Sigma[2,2] = sqrt(1.0 - rho_w^2);
  
  matrix[N1, N1] L_K = build_cholesky_factor(W1[,1], loc[1:N1,], 1.0, rho, delta);
}

model{
  
  vector[2*N1] W1_vec = to_vector(W1);
  vector[2*N1] mu_vec = to_vector(mu_mat);
  
    mu_x ~ normal(0, 5);
    sigma_x ~ inv_gamma(2,1);
    rho_x ~ inv_gamma(2,1);
    alpha_x ~ normal(0,5);

    for (p in 1:P) {
      matrix[N1, N1] L_K_x = build_cholesky_factor(X1[,p], loc[1:N1,], alpha_x[p], rho_x[p], sigma_x[p]);
      vector[N1] mu_p = rep_vector(mu_x[p], N1);
      
      X1[,p] ~ multi_normal_cholesky(mu_p, L_K_x);
    }
    
    
    mu0 ~ normal(0, 10);
    to_vector(coef_mat) ~ normal(0, 10);
    rho ~ inv_gamma(2,1);
    sigma_w_sq ~ inv_gamma(2,1);
    rho_w ~ std_normal();
    // to_vector(Winc) ~ normal(0, 10);
    to_vector(eta) ~ std_normal();
    xi ~ std_normal();
    latent_lengths ~ normal(0, 10);
    
    matrix[2,2] Sigma_inv;
    Sigma_inv = chol2inv(L_Sigma);
  
    matrix[N1, N1] K_inv = chol2inv(L_K);
  
    // matrix[2*N1, 2*N1] Sigma_inv_K_inv_kprod;
    // Sigma_inv_K_inv_kprod = kronecker_prod(Sigma_inv, K_inv);
    
    // mu fc
  // matrix[2*P,2*P] Q_beta = X_aug' * Sigma_inv_K_inv_kprod * X_aug + 1.0/lambda0^2 * identity_matrix(2*P);
  // matrix[2*P,2*P] L_Q = cholesky_decompose(Q_beta);
  // matrix[2*P,2*P] L_Q_inv = mdivide_left_tri_low(L_Q, identity_matrix(2*P));
  // 
  // vector[2*P] E_beta;
  // E_beta = chol2inv(L_Q) * X_aug' * Sigma_inv_K_inv_kprod * W1_vec;
  // 
  // to_vector(B) ~ multi_normal_cholesky(E_beta, L_Q_inv');
  // 
  // lengths fc
  vector[N_obs] lengths_llpd = latent_length_llpd(latent_lengths, Wobs, mu_mat[1:N_obs,], K_inv[1:N_obs, 1:N_obs], Sigma_inv, N_obs);
  target += sum(lengths_llpd);
  
  // // cov fc + spatial Par fc
  // target += (-N1 -2*a_s - 2)* log(sigma_w) - N1/2.0 * log(1 - rho_w^2) - sum(log(diagonal(L_K))) - (W1_vec - mu_vec)' * Sigma_inv_K_inv_kprod * (W1_vec - mu_vec) - b_s / sigma_w^2;
    
    matrix[N1, 2] f = L_K * eta * diag_pre_multiply(alpha, L_Sigma)';

    to_vector(W1) ~ normal(to_vector(mu_mat + f), 1.0);
}

generated quantities{
    real spat_decay = 1.0 / rho;
    
    matrix[P, 2] B;
    B = R_ast_inverse * coef_mat; // coefficients on x
    
    // In-Sample PPD Draws
    matrix[N1, P] x_pred1;
      
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

    matrix[N1, 2] y_pred1;
    vector[N1] theta_pred1;
    {
      matrix[N1, 2] f_1;
      // matrix[N1, N1] K = gp_exponential_cov(loc[1:N1,], 1.0, rho);
      // matrix[N1, N1] L_K;
      // 
      // // diagonal elements
      // for (n in 1:N1) {
      //   K[n, n] = K[n, n] + delta;
      // }
      // 
      // L_K = cholesky_decompose(K);
      f_1 = L_K * eta * diag_pre_multiply(alpha, L_Sigma)';
    
      for(i in 1:N1) {
          // y_pred1[i,] = multi_normal_rng(mu_mat[i,]' + f_1[i,]', identity_matrix(2))';
          y_pred1[i,1] = std_normal_rng() + mu_mat[i,1] + f_1[i,1];
          y_pred1[i,2] = std_normal_rng() + mu_mat[i,2] + f_1[i,2];
          theta_pred1[i] = atan2(y_pred1[i,2], y_pred1[i,1]);
          // theta_pred1[i] = atan2(W1[i,2], W1[i,1]);
      }
    }
    
    vector[N_inc] theta_imp;
    for (n in 1:N_inc) {
      theta_imp[n] = atan2(Winc[n,2], Winc[n,1]);
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

    matrix[N2, 2] y_pred2;
    vector[N2] theta_pred2;
    {
      vector[2*N2] vec_Y2;
      // matrix[N1, N1] K = gp_exponential_cov(loc[1:N1,], 1.0, rho);
      // matrix[N1, N1] L_K;
      // // diagonal elements
      // for (n in 1:N1) {
      //   K[n, n] = K[n, n] + delta;
      // }
      // L_K = cholesky_decompose(K);

      matrix[N2, N2] K_pred = gp_exponential_cov(loc2, 1.0, rho);
      K_pred = add_diag(K_pred, delta);
      
      matrix[N2, N2] L_Kpred = cholesky_decompose(K_pred);

      matrix[N1, N2] K_loc1_loc2 = gp_exponential_cov(loc[1:N1,], loc2, 1.0, rho);

      vector[2*N1] v1_minus_mu = to_vector(W1' - mu_mat');
      matrix[N1, N2] U = mdivide_left_tri_low(L_K, K_loc1_loc2);

      matrix[N2, 2] mu2 = rep_matrix(mu0', N2) + x_pred2 * B;
      vector[2*N2] cond_mean;
      matrix[2*N2, 2*N2] cond_cov;

      cond_mean = to_vector(mu2') + kronecker_prod(mdivide_right_tri_low(U', L_K), identity_matrix(2)) * v1_minus_mu;

      matrix[2,2] Sigma = tcrossprod(L_Sigma);
      cond_cov = kronecker_prod(K_pred - crossprod(U), Sigma);

      vec_Y2 = cond_mean + cholesky_decompose(cond_cov) * xi;
      y_pred2 = to_matrix(vec_Y2, 2, N2)';
    }


    for(i in 1:N2) {
        theta_pred2[i] = atan2(y_pred2[i,2], y_pred2[i,1]);
    }
}


