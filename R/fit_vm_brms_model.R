#' Fit the von Mises Regression Using Stan
#' @export
#' @family models
#' @param theta vector of angles
#' @param X vector or matrix of predictors
#' @param X_ppd vector or matrix of holdout predictors
#' @param iter_sampling int number of sampling iterations
#' @param iter_warmup int number of warmup iterations
#' @param refresh int how often to print updates to console
#' @param ... 
#'
#' @return CmdStan fit object
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
#'   fit_vm_brms_model(theta = df$theta, X = df$X, X_ppd = df$X)
#' }
fit_vm_brms_model <- function(theta, X, X_ppd, iter_sampling = 1000, iter_warmup = 1000, refresh = 500,...) {
  stopifnot(is.numeric(theta) && max(abs(theta)) < 2*pi)
  
  X_mat <- stats::model.matrix(~., data = as.data.frame(X))
  X_ppd_mat <- stats::model.matrix(~., data = as.data.frame(X_ppd))
  
  model <- instantiate::stan_package_model(
    name = "vm_brms_reg",
    package = "pnregstan"
  )
  data_list <- list(N = length(theta),
                    N_ppd = nrow(X_ppd_mat),
                    K = ncol(X_mat),
                    Y = theta,
                    X = X_mat,
                    X_ppd = X_ppd_mat, 
                    prior_only = 0)
  fit <- model$sample(data = data_list, 
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = 2,
                      ...)
  
  return(fit)
}
