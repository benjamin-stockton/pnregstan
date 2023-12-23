#' Simulate Projected Normal Data from a latent ARX process
#'
#' @param N int time series length
#' @param mu_0  initial latent mean vector
#' @param B_vec vector of regression coefficients
#' @param Psi_vec vector of autoregressive coefficients
#' @param Sigma_vec vector of stationary covariance parameters
#' @param X NULL or matrix of predictor time series
#' @param ar_X1 NULL or vector of autoregressive coefficients for simulated X1
#' @param ma_X1 NULL or vector of moving average coefficients for simulated X1 
#' @param ar_X2 NULL or vector of autoregressive coefficients for simulated X2
#' @param ma_X2 NULL or vector of moving average coefficients for simulated X2 
#'
#' @return df a data frame with theta, U1, U2, and either simulated X time series or provided X_mat time series
#' @export
#'
#' @examples
#' df <- pn_arx_sim_data(N = 100, ar_X1 = c(0.24), ma_X1 = numeric(0),
#'                       ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
pn_arx_sim_data <- function(N = 100, mu_0 = c(0,0), B_vec = c(1, .3, 0, -2.5), Psi_vec = c(0.75, -.2, -.2, .5), Sigma_vec = c(1,0,0,1), X = NULL, ar_X1 = NULL, ma_X1 = NULL, ar_X2 = NULL, ma_X2 = NULL) {
  
  B <- matrix(B_vec, nrow = 2)
  Psi <- matrix(Psi_vec, nrow = 2)
  Sigma <- matrix(Sigma_vec, nrow = 2)
  
  if (is.null(ar_X1) && is.na(ar_X1) && is.null(X)) {
    ar_X1 <- c(0.24)
    ma_X1 <- numeric(0)
  }
  if (is.null(ar_X2) && is.na(ar_X2) && is.null(X)) {
    ar_X2 <- c(0.75)
    ma_X2 <- c(-0.25, 0.25)
  }
  
  if (is.null(X)) {
    X_mat <- cbind(
      stats::arima.sim(n = N, list(ar = ar_X1, ma = ma_X1)),
      stats::arima.sim(n = N, list(ar = ar_X2, ma = ma_X2))
    )
  }
  else {
    X_mat <- X
  }
  
  mu <- X_mat %*% B
  
  Y <- matrix(NA, nrow = N, ncol = 2)
  U <- matrix(NA, nrow = N, ncol = 2)
  eps_0 <- mvtnorm::rmvnorm(1, rep(0,2), Sigma)
  Y[1,] <- mu_0 + mu[1,] + eps_0
  U[1,] <- Y[1,] / sqrt(Y[1,1]^2 + Y[1,2]^2)
  for (i in 2:N) {
    eps_i <- mvtnorm::rmvnorm(1, rep(0,2), Sigma)
    Y[i,] <- mu_0 + Y[i-1,] %*% Psi + mu[i,] + eps_i
    U[i,] <- Y[i,] / sqrt(Y[i,1]^2 + Y[i,2]^2)
  }
  
  theta <- atan2(U[,2], U[,1])
  
  df <- data.frame(
    theta = theta,
    U1 = U[,1],
    U2 = U[,2]
  )
  
  df <- cbind(df, X_mat)
  if (is.data.frame(X_mat)) {
    colnames(df) <- c("theta", "U1", "U2", colnames(X_mat))
  }
  else {
    colnames(df) <- c("theta", "U1", "U2", paste0("X", 1:ncol(X_mat)))
  }
  return(df)
}