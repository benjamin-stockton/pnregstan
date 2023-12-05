#' Simulate Projected Normal Regression Data
#'
#' @param N Sample size
#' @param B Regression coefficients passed as a vector turned into a matrix (Px2) by row
#' @param Sigma_theta Projected normal covariance values passed as a vector turned into a matrix (2x2) by row, should be symmetric
#' @param X a matrix (or data frame) of (N x P) predictors or NULL
#' @param mu_X numeric for mean of the distribution of a single simulated predictor X or NULL
#' @param sigma_X non-negative numeric for the standard deviation of the single simulated predictor X or NULL
#'
#' @return df A data frame with three columns for theta, U1, U2, and columns for X (1 col simulated or the passed X appended)
#' @export
#'
#' @examples
#'  df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
pnreg_sim_data <- function(N = 100, B = c(5, 0.5, 0.75, 0.25), Sigma_theta = c(1,0,0,1), X = NULL, mu_X = NULL, sigma_X = NULL) {
  B <- matrix(B, ncol = 2, byrow = TRUE)
  Sigma_theta <- matrix(Sigma_theta, nrow = 2, ncol = 2, byrow = TRUE)
  if (is.null(X) && is.numeric(mu_X) && is.numeric(sigma_X)) {
    X <- stats::rnorm(N, mu_X, sigma_X)
    
    X_mat <- stats::model.matrix(~., data = data.frame(X = X))
  }
  else {
    X_mat <- stats::model.matrix(~., data = as.data.frame(X))
  }
  
  mu <- X_mat %*% B
  U <- matrix(NA, nrow = N, ncol = 2)
  
  for (i in 1:N) {
    tmp <- mvtnorm::rmvnorm(1, mu[i,], Sigma_theta)
    U[i,] <- tmp / sqrt(sum(tmp^2))
  }
  
  theta <- atan2(U[,2], U[,1])
  
  df <- data.frame(
    theta = theta,
    U1 = U[,1],
    U2 = U[,2]
  )
  df <- cbind(df, X)
  
  return(df)
}