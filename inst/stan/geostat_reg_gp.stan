functions {
    matrix build_cholesky_factor(vector Y1,
                                 array[] vector loc1,
                                 real alpha,
                                 real rho,
                                 real sigma) {
        int N1 = rows(Y1);
        
        matrix[N1, N1] K = gp_exp_quad_cov(loc1, alpha, rho);
        real sq_sigma = square(sigma);
        
        for (i in 1:N1) {
            K[i, i] = K[i, i] + sq_sigma;
        }
        
        matrix[N1, N1] L_K = cholesky_decompose(K);
        return L_K;
    }
    
    vector gp_pred_rng(array[] vector loc2,
                     vector Y1,
                     array[] vector loc1,
                     vector mu1,
                     vector mu2,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = rows(Y1);
    int N2 = size(loc2);
    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] y_minus_mu;
      vector[N1] K_div_Y1;
      matrix[N1, N2] k_X1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      K = gp_exp_quad_cov(loc1, alpha, rho);
      for (n in 1:N1) {
        K[n, n] = K[n, n] + square(sigma);
      }
      L_K = cholesky_decompose(K);
      y_minus_mu = Y1 - mu1;
      K_div_Y1 = mdivide_left_tri_low(L_K, y_minus_mu);
      K_div_Y1 = mdivide_right_tri_low(K_div_Y1', L_K)';
      k_X1_x2 = gp_exp_quad_cov(loc1, loc2, alpha, rho);
      f2_mu = mu2 + (k_X1_x2' * K_div_Y1);
      v_pred = mdivide_left_tri_low(L_K, k_X1_x2);
      cov_f2 = gp_exp_quad_cov(loc2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
 int<lower=1> N1; //number of data points
 int<lower=1> N2;
 vector[N1] Y1; //number of outcomes
 array[N1] vector[2] loc1;
 array[N2] vector[2] loc2;
 int P; // number of predictors
 matrix[N1, P] X1; // predictor
 // int<lower=1> N2; // number of new points
 // matrix[N1+N2, N1+N2] dist; //distances between points
}

transformed data{
 // int<lower=1> N;
 // N = N1 + N2;
 vector[N1] mu_0 = rep_vector(0, N1);
 real delta = 1e-9;
}

parameters{
 // vector[N1] z1;
 // vector[N2] z2;
 // array[N2] vector[P] z2;
 vector[P] beta;
 real beta0;
 real<lower=0> sigma;
 real<lower=0>  rho;
 real<lower=0> alpha;
 
 // vector[P] mu_x;
 vector<lower=0>[P] sigma_x;
 vector<lower=0>[P]  rho_x;
 vector<lower=0>[P] alpha_x;
}

transformed parameters{
 vector[N1] mu;
 mu = rep_vector(beta0, N1);
 for(i in 1:N1) {
     mu[i] += X1[i,] * beta;
 }
}

model{
    // mu_x ~ normal(0, 100);
    sigma_x ~ normal(0,5);
    rho_x ~ inv_gamma(5,5);
    alpha_x ~ normal(0,5);
      
    for (p in 1:P) {
      matrix[N1, N1] L_K = build_cholesky_factor(X1[,p], loc1, alpha_x[p], rho_x[p], sigma_x[p]);
      
      X1[,p] ~ multi_normal_cholesky(mu_0, L_K);
    }
    
  
    sigma ~ normal(0, 5);
    rho ~ inv_gamma(5,5);
    beta0 ~ normal(0, 5);
    beta ~ normal(0,5);
    alpha ~ std_normal();
  
    matrix[N1,N1] L = build_cholesky_factor(Y1, loc1, alpha, rho, sigma);
 
    Y1 ~ multi_normal_cholesky(mu, L);
}

generated quantities{
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
    
    vector[N1] y_pred1;
    vector[N1] f1;
    f1 = gp_pred_rng(loc1, 
                     Y1,
                     loc1,
                     mu,
                     mu,
                     alpha,
                     rho,
                     sigma,
                     delta);
    for(i in 1:N1) {
        y_pred1[i] = normal_rng(f1[i], sigma);
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
    
    vector[N2] f2;
    vector[N2] y_pred2;
    vector[N2] mu_pred2 = rep_vector(beta0, N2) + x_pred2 * beta;
    f2 = gp_pred_rng(loc2, 
                     Y1,
                     loc1,
                     mu,
                     mu_pred2,
                     alpha,
                     rho,
                     sigma,
                     delta);
    for (i in 1:N2) {
        y_pred2[i] = normal_rng(f2[i], sigma);
    }
}
