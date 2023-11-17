#' @title Fit the Projected Normal Regression model with Identity Covariance.
#' @export
#' @family models
#' @description Fit the PN Regression Stan model and return posterior summaries.
#' @return A data frame of posterior summaries.
#' @param theta a vector length N of angles
#' @param X a matrix of size N x P of predictors
#' @param X_ppd a matrix of size N_tilde x P of predictors for the posterior predictive draws
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   fit_pnreg_identity_model(theta, X, X_ppd)
#' }
fit_pnreg_identity_model <- function(theta, X, X_ppd, ...) {
  stopifnot(is.numeric(theta) && range(theta) < 2*pi)
  U <- cbind(cos(theta), sin(theta))
  
  model <- instantiate::stan_package_model(
    name = "pnreg_identity",
    package = "pnregstan"
  )
  data_list <- list(N = length(theta),
                    N_tilde = nrow(X_ppd),
                    P = ncol(X),
                    U = U,
                    X = X,
                    X_tilde = X_ppd)
  fit <- model$sample(data = data_list, ...)
  
  return(fit)
}
