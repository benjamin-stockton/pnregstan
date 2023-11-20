#' @title Fit the Projected Normal Regression model with General Covariance.
#' @export
#' @family models
#' @description Fit the PN Regression Stan model and return posterior summaries. Covariance follows the constrained parameterization from 
#' @return A data frame of posterior summaries.
#'
#' @param theta a vector length N of angles
#' @param X a matrix of size N x P of predictors
#' @param iter_sampling int number of iterations to sample per chain
#' @param iter_warmup int number of warmup iterations to sample per chain
#' @param refresh int how often CmdStan should print an update. Set to 0 to suppress updates.
#' @param X_ppd a matrix of size N_tilde x P of predictors for the posterior predictive draws
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
#'   fit_pnreg_identity_model(theta = df$theta, X = df$X, X_ppd = df$X)
#' }
fit_pnreg_gen_model <- function(theta, X, X_ppd, iter_sampling = 1000, iter_warmup = 3000, refresh = 500) {
  stopifnot(is.numeric(theta) && max(abs(theta)) < 2*pi)
  U <- angle_to_unit_vec(theta)
  X_mat <- stats::model.matrix(~., data = as.data.frame(X))
  X_ppd_mat <- stats::model.matrix(~., data = as.data.frame(X_ppd))
  
  model <- instantiate::stan_package_model(
    name = "proj_normal_general_circ_reg",
    package = "pnregstan"
  )
  data_list <- list(N = length(theta),
                    N_ppd = nrow(X_ppd_mat),
                    P = ncol(X_mat),
                    U = U,
                    X = X_mat,
                    X_ppd = X_ppd_mat)
  fit <- model$sample(data = data_list, 
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      show_exceptions = FALSE,
                      show_messages = FALSE,
                      adapt_delta = 0.85,
                      chains = 2)
  
  return(fit)
}
