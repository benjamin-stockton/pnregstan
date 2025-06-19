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
 int<lower=0> P;
 vector[N1] theta1; // angular outcomes
 matrix[N1, P] X1;
 array[N1] vector[2] loc1;
 array[N2] vector[2] loc2;
}

transformed data{
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
 vector[2] alpha = rep_vector(1.0, 2);
 
 vector[N1] zeros = rep_vector(0, N1);

 // matrix[2*N1, 2*P] X_aug = kronecker_prod(identity_matrix(2), X1);
 
 // hyperparameters
 // real<lower=0> lambda0 = 100;
 // real<lower=0> a_s = 2;
 // real<lower=0> b_s = 1;
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
 
  // vector[P] mu_x;
 vector<lower=0>[P] sigma_x;
 vector<lower=0>[P]  rho_x;
 vector<lower=0>[P] alpha_x;
}

transformed parameters {
  matrix[N1, 2] W1;
  
  W1 = diag_pre_multiply(latent_lengths, U1);
  
  matrix[N1, 2] mu_mat;
  for (i in 1:N1) {
    mu_mat[i,] = mu0' + X1[i,] * B;
  }
  
  matrix[2,2] L_Sigma;
  L_Sigma[1,1] = sigma_w;
  L_Sigma[1,2] = 0;
  L_Sigma[2,1] = rho_w;
  L_Sigma[2,2] = sqrt(1.0 - rho_w^2);
}

model{
  
    sigma_x ~ normal(0,5);
    rho_x ~ inv_gamma(2,1);
    alpha_x ~ normal(0,5);
      
    for (p in 1:P) {
      matrix[N1, N1] L_K = build_cholesky_factor(X1[,p], loc1, alpha_x[p], rho_x[p], sigma_x[p]);
      
      X1[,p] ~ multi_normal_cholesky(zeros, L_K);
    }
    
    
    mu0 ~ normal(0, 10);
    to_vector(B) ~ normal(0, 10);
    rho ~ inv_gamma(1.1, 1);
    sigma_w ~ inv_gamma(1.1, 1);
    rho_w ~ std_normal();
    // alpha ~ std_normal();
    // L_Sigma ~ lkj_corr_cholesky(3);
    to_vector(eta) ~ std_normal();
    latent_lengths ~ normal(0, 10);
    
    matrix[N1, 2] f;
    matrix[2,2] Sigma_inv;
    Sigma_inv = chol2inv(L_Sigma);
  
    matrix[N1, N1] L_K = build_cholesky_factor(W1[,1], loc1, 1.0, rho, delta);
    matrix[N1, N1] K_inv = chol2inv(L_K);
    {
      f = L_K * eta * diag_pre_multiply(alpha, L_Sigma)';
    }

    // lengths fc
    vector[N1] lengths_llpd = latent_length_llpd(latent_lengths, W1, mu_mat, K_inv, Sigma_inv, N1);
    target += sum(lengths_llpd);

    to_vector(W1) ~ normal(to_vector(mu_mat + f), 1.0);
}

generated quantities{
    real spat_decay = 1.0/rho;
  
    // In-Sample PPD Draws
    matrix[N1, P] x_pred1;
    vector[N1] mu_pred_0 = rep_vector(0, N1);
    for (p in 1:P) {
      
      vector[N1] f1_x;
      f1_x = gp_pred_rng(loc1, 
                     X1[,p],
                     loc1,
                     mu_pred_0,
                     mu_pred_0,
                     alpha_x[p],
                     rho_x[p],
                     sigma_x[p],
                     delta);
      
      x_pred1[,p] = multi_normal_rng(f1_x, square(sigma_x[p]) * identity_matrix(N1));
      
    }

    matrix[N1, 2] y_pred1;
    vector[N1] theta_pred1;
    {
      matrix[N1, 2] f_1;
      matrix[N1, N1] K = gp_exponential_cov(loc1, 1.0, rho);
      matrix[N1, N1] L_K;

      // diagonal elements
      for (n in 1:N1) {
        K[n, n] = K[n, n] + delta;
      }

      L_K = cholesky_decompose(K);
      f_1 = L_K * eta * diag_pre_multiply(alpha, L_Sigma)';
    
      for(i in 1:N1) {
          y_pred1[i,] = multi_normal_rng(mu_mat[i,]' + f_1[i,]', identity_matrix(2))';
          theta_pred1[i] = atan2(y_pred1[i,2], y_pred1[i,1]);
      }
    }
    
    // Regular Grid/Out-of-Sample PPD Draws
    matrix[N2,P] x_pred2;
    vector[N2] mu_pred_02 = rep_vector(0, N2);
    for (p in 1:P) {
      vector[N2] f2_x;
      f2_x = gp_pred_rng(loc2, 
                     X1[,p],
                     loc1,
                     mu_pred_0,
                     mu_pred_02,
                     alpha_x[p],
                     rho_x[p],
                     sigma_x[p],
                     delta);
                     
      x_pred2[,p] = multi_normal_rng(f2_x, square(sigma_x[p]) * identity_matrix(N2));
    }

    matrix[N2, 2] y_pred2;
    vector[N2] theta_pred2;
    {
      vector[2*N2] vec_Y2;
      vector[2*N2] xi;
      matrix[N1, N1] K = gp_exponential_cov(loc1, 1.0, rho);
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
      
      matrix[N1, N2] K_loc1_loc2 = gp_exponential_cov(loc1, loc2, 1.0, rho);
      
      vector[2*N1] v1_minus_mu = to_vector(W1' - mu_mat');
      matrix[N1, N2] U = mdivide_left_tri_low(L_K, K_loc1_loc2);
      
      matrix[N2, 2] mu2 = rep_matrix(mu0', N2) + x_pred2 * B;
      vector[2*N2] cond_mean;
      matrix[2*N2, 2*N2] cond_cov;
      
      cond_mean = to_vector(mu2') + kronecker_prod(mdivide_right_tri_low(U', L_K), identity_matrix(2)) * v1_minus_mu;
      
      matrix[2,2] Sigma = L_Sigma*L_Sigma';
      cond_cov = kronecker_prod(K_pred - crossprod(U), Sigma);
      
      vec_Y2 = cond_mean + cholesky_decompose(cond_cov) * xi;
      y_pred2 = to_matrix(vec_Y2, 2, N2)';
    }
    
    
    for(i in 1:N2) {
        theta_pred2[i] = atan2(y_pred2[i,2], y_pred2[i,1]);
    }
}


