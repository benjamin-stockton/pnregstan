#' Projected Normal AR Process with Exogenous Predictors (Identity Covariance)
#' @description
#' The model is fitted using the projected normal regression with identity covariance.
#' 
#' @param theta vector angles
#' @param X matrix of exogenous predictors
#' @param X_ppd matrix of exogenous predictors for PPD draws
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return vector of ppd drawn angles
#' @export
#'
#' @examples
#' ## Not run:
#' if (instantiate::stan_cmdstan_exists()) {
#'   df <- pn_arx_sim_data(N = 100, 
#'   ar_X1 = c(0.24), ma_X1 = numeric(0), 
#'   ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
#'   
#'   df_ppd <- df[91:100,]
#'   fit_pn_arx_model(theta = df$theta[1:90], 
#'   X = df[1:90,c("X1", "X2")], 
#'   X_ppd = df_ppd[,c("X1", "X2")])
#' }
#' ## End(Not run)
fit_pn_arx_id_model <- function(theta, X, X_ppd, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  stopifnot(is.numeric(theta) && max(abs(theta)) < 2*pi)
  U <- angle_to_unit_vec(theta)
  
  model <- instantiate::stan_package_model(
    name = "pn_arx_id_process",
    package = "pnregstan"
  )
  data_list <- list(N = length(theta),
                    N_ppd = nrow(X_ppd),
                    K = ncol(X),
                    U = U,
                    X = X,
                    X_ppd = X_ppd, 
                    sigma_0 = 10)
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)

  return(fit)
}