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
  int<lower=1> N1; //number of data points
  int<lower=1> N2;
  vector[N1] theta1; // angular outcomes
  vector[N1] mag1; // magnitudes for angles
  array[N1] vector[2] loc1;
  array[N2] vector[2] loc2;
}

transformed data {
 int<lower=1> N;
 N = N1 + N2;
 real delta = 1e-9;
 matrix[N1, 2] U1;
 for (i in 1:N1) {
   U1[,1] = cos(theta1);
   U1[,2] = sin(theta1);
 }
 
 array[N] vector[2] loc;
 for (n in 1:N1) {
   loc[n] = loc1[n];
 }
 for (n in 1:N2) {
   loc[N1+n] = loc2[n];
 }
 // vector[2] mu = rep_vector(0, 2);
 // matrix[2,2] L_Sigma = identity_matrix(2);
 // vector[2] alpha = rep_vector(1.0, 2);
 
 matrix[N1,1] ones = to_matrix(rep_vector(1.0, N1), N1, 1);
 matrix[2*N1, 2] A = kronecker_prod(identity_matrix(2), ones);
 
 // hyperparameters
 real<lower=0> lambda0 = 100;
 real<lower=0> a_s = 2.0;
 real<lower=0> b_s = 1.0;
 real<lower=0> a_ma = 2.0;
 real<lower=0> b_ma = 1.0;
 real<lower=0> a_ms = 2.0;
 real<lower=0> b_ms = 1.0;
 
 matrix[N1, N1] dist_mat;
 dist_mat = -log(gp_exponential_cov(loc1, 1.0, 1.0));
 real max_dist = max(to_vector(dist_mat));
 
 // real rho_w = 0.0;
 // real sigma_w = 1.0;
 // real rho = 0.5;
 // real alpha_m = 1.0;
}

parameters {
 vector<lower=0>[N1] latent_lengths;
 vector[2] mu_w;
 
 real<lower=0> rho;
 
 real<lower=0> sigma_w;
 real<lower=-1, upper=1> rho_w;
 
 vector[3] beta_m;
 
 real<lower=0> sigma_m;
 real<lower=0> alpha_m;
 real<lower=0> rho_m;
}

transformed parameters {
  real sigma_w_sq = sigma_w^2;
  real sigma_m_sq = sigma_m^2;
  real phi = 1.0/rho;
  real phi_m = 1.0/rho_m;
  matrix[N1, 2] W1;
  W1 = diag_pre_multiply(latent_lengths, U1);
  
  matrix[2,2] L_Sigma;
  L_Sigma[1,1] = sigma_w;
  L_Sigma[1,2] = 0;
  L_Sigma[2,1] = rho_w;
  L_Sigma[2,2] = sqrt(1.0 - rho_w^2);
}

model {
  latent_lengths ~ normal(0, 10);
  
  vector[2*N1] W1_vec = to_vector(W1);
  
  matrix[2,2] Sigma_inv;
  Sigma_inv = chol2inv(L_Sigma);
  
  matrix[N1, N1] L_K = build_cholesky_factor(theta1, loc1, 1.0, rho, delta);
  matrix[N1, N1] K_inv = chol2inv(L_K);
  
  matrix[2*N1, 2*N1] Sigma_inv_K_inv_kprod;
  Sigma_inv_K_inv_kprod = kronecker_prod(Sigma_inv, K_inv);
  
  vector[N1] beta_4_s;
  beta_4_s = beta_m[2] * U1[,1] + beta_m[3] * U1[,2];
  
  matrix[N1, N1] L_K_m = build_cholesky_factor(mag1, loc1, alpha_m, rho_m, sigma_m);
  
  matrix[N1, N1] Q_mag = chol2inv(L_K_m);
  matrix[N1, N1] K_m = tcrossprod(L_K_m);
  
  // Priors
  mu_w ~ normal(0, lambda0);
  sigma_w_sq ~ inv_gamma(a_s, b_s);
  // sigma_w_sq ~ normal(0, 10);
  rho_w ~ uniform(-1,1);
  // rho ~ uniform(0, max_dist);
  phi ~ uniform(0, max_dist);
  // rho ~ normal(0, max_dist / 3);
  
  // mu fc
  // matrix[2,2] Q_mu = A' * Sigma_inv_K_inv_kprod * A + 1.0/lambda0^2 * identity_matrix(2);
  // matrix[2,2] L_Q = cholesky_decompose(Q_mu);
  // matrix[2,2] L_Q_inv = mdivide_left_tri_low(L_Q, identity_matrix(2));
  // 
  // vector[2] E_mu;
  // E_mu = chol2inv(L_Q_inv) * A' * Sigma_inv_K_inv_kprod * W1_vec;
  // 
  // mu_w ~ multi_normal_cholesky(E_mu, L_Q_inv');
  
  // lengths fc
  vector[N1] lengths_llpd;
  {
    real a;
    real b;
    vector[2] E_si;
    matrix[2,2] w_ii;
    vector[N1] mag1_c = mag1 - beta_m[1];
    for (i in 1:N1) {

      w_ii = K_inv[i,i] * tcrossprod(L_Sigma);
      E_si = mu_w - w_ii * Sigma_inv * (W1 - rep_matrix(mu_w', N1))' * K_inv[,i] + w_ii * K_inv[i,i] * (W1[i,] - mu_w')';
      real E_si2 = mag1_c[i] + (K_m[i,i]) * sum(Q_mag[,i] .* (mag1_c - latent_lengths .* beta_4_s)) - (mag1[i] - beta_m[1] - latent_lengths[i] * beta_4_s[i]);

      a = quad_form(Sigma_inv, U1[i,]');
      b = U1[i,] * Sigma_inv * E_si;
      lengths_llpd[i] = log(latent_lengths[i]) - K_inv[i,i] / 2.0 * (a + beta_4_s[i]^2) * (latent_lengths[i] - (b + beta_4_s[i] * E_si2) / (a + beta_4_s[i]^2))^2;
      // lengths_llpd[i] = log(latent_lengths[i]) - K_inv[i,i] / 2.0 * (latent_lengths[i]^2 * a - 2 * latent_lengths[i] * b);

    }
  }

  target += sum(lengths_llpd);
  
  // cov fc + spatial Par fc
  target += (-N1 -2*a_s - 2)* log(sigma_w) - N1/2.0 * log(1 - rho_w^2) - sum(log(diagonal(L_K))) - (W1_vec - A * mu_w)' * Sigma_inv_K_inv_kprod * (W1_vec - A * mu_w) - b_s / sigma_w^2;
  
  // W1_vec ~ multi_normal_prec(A * mu_w, Sigma_inv_K_inv_kprod);
  
  // Magnitude Model
  
  beta_m ~ normal(0, 10);
  sigma_m_sq ~ inv_gamma(a_ms,b_ms);
  // sigma_m ~ normal(0, 10);
  // phi_m ~ uniform(0, max_dist);
  rho_m ~ uniform(0, max_dist);
  // rho_m ~ normal(0, max_dist/3.0);
  alpha_m ~ inv_gamma(a_ma,b_ma);
  // alpha_m ~ normal(0, 10);
  
  vector[N1] mean_m = beta_m[1] + W1 * beta_m[2:3];
  // vector[N1] res_m = mag1 - mean_m;
  // 
  // matrix[N1, 3] W1_ones = append_col(rep_vector(1.0, N1), W1);
  // matrix[3, N1] L_beta = mdivide_right_tri_low(W1_ones', L_K_m);
  // matrix[3, 3] Q_beta = tcrossprod(L_beta);
  // vector[3] E_beta = inverse(Q_beta) * W1_ones' * Q_mag * mag1;
  // 
  // beta_m ~ multi_normal_prec(E_beta, Q_beta);
  // 
  // target += -log(determinant(Q_mag)) - quad_form(Q_mag, res_m) - 2*(a_ma + 1) * log(alpha_m) - b_ma / alpha_m^2 - 2*(a_ms + 1) * log(sigma_m) - b_ms / sigma_m^2;
  
  mag1 ~ multi_normal_cholesky(mean_m, L_K_m);
}

generated quantities{
    real spat_decay = 1.0/rho;
    real mean_dir = atan2(mu_w[2], mu_w[1]);
  
    // In-Sample PPD Draws
    matrix[N1, 2] y_pred1;
    vector[N1] theta_pred1;
    
    vector[N1] mag_pred1;
    vector[N1] mu_pred_m;
    mu_pred_m = beta_m[1] + W1 * beta_m[2:3];
    matrix[N1, N1] L_K = build_cholesky_factor(theta1, loc1, 1.0, rho, delta);
    {
      matrix[N1, 2] f_1;
      matrix[N1, 2] eta;
      
      for (n in 1:N1) {
        eta[n,1] = normal_rng(0,1);
        eta[n,2] = normal_rng(0,1);
      }
      
      f_1 = L_K * eta * L_Sigma';
      
      // GP for the magnitudes
      vector[N1] f1_m;
      f1_m = gp_pred_rng(loc1, 
                     mag1,
                     loc1,
                     mu_pred_m,
                     mu_pred_m,
                     alpha_m,
                     rho_m,
                     sigma_m,
                     delta);
      for (i in 1:N1) {
        y_pred1[i,] = mu_w' + f_1[i,];
        theta_pred1[i] = atan2(y_pred1[i,2], y_pred1[i,1]);
        mag_pred1[i] = normal_rng(f1_m[i], sigma_m);
      }
      
    }
    
    // Regular Grid/Out-of-Sample PPD Draws
    
    matrix[N2, 2] y_pred2;
    vector[N2] theta_pred2;
    {
      vector[2*N2] vec_Y2;
      vector[2*N2] xi;
      
      matrix[N2, N2] K_pred = gp_exponential_cov(loc2, 1.0, rho);
      
      K_pred = add_diag(K_pred, delta);
      for (n in 1:N2) {
        // K_pred[n,n] += delta;
        
        // draws from standard normal
        xi[n] = std_normal_rng();
        xi[N2+n] = std_normal_rng();
      }
      matrix[N2, N2] L_Kpred = cholesky_decompose(K_pred);
      
      matrix[N1, N2] K_loc1_loc2 = gp_exponential_cov(loc1, loc2, 1.0, rho);
      
      vector[2*N1] v1_minus_mu = to_vector(W1' - rep_matrix(mu_w', N1)');
      matrix[N1, N2] U = mdivide_left_tri_low(L_K, K_loc1_loc2);
      
      matrix[N2, 2] mu2 = rep_matrix(mu_w', N2);
      vector[2*N2] cond_mean;
      matrix[2*N2, 2*N2] cond_cov;
      
      cond_mean = to_vector(mu2') + kronecker_prod(mdivide_right_tri_low(U', L_K), identity_matrix(2)) * v1_minus_mu;
      
      matrix[2,2] Sigma = L_Sigma*L_Sigma';
      cond_cov = kronecker_prod(K_pred - crossprod(U), Sigma);
      
      vec_Y2 = cond_mean + cholesky_decompose(cond_cov) * xi;
      y_pred2 = to_matrix(vec_Y2, 2, N2)';
    }
    
    vector[N2] mag_pred2;
    vector[N2] mu_pred_m2;
    mu_pred_m2 = beta_m[1] + y_pred2 * beta_m[2:3];
    {
      vector[N2] f2_m;
      f2_m = gp_pred_rng(loc2, 
                     mag1,
                     loc1,
                     mu_pred_m,
                     mu_pred_m2,
                     alpha_m,
                     rho_m,
                     sigma_m,
                     delta);
                     
      for (i in 1:N2) {
        mag_pred2[i] = normal_rng(f2_m[i], sigma_m);
        theta_pred2[i] = atan2(y_pred2[i,2], y_pred2[i,1]);
      }
    }
}

