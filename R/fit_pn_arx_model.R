#' Title
#'
#' @param theta vector angular time series; values should be between 0 and 2*pi
#' @param X vector or matrix of predictor time series for training data
#' @param X_ppd vector or matrix of predictor time series for testing data
#' @param iter_sampling int MCMC sampling iterations to send to CmdStan
#' @param iter_warmup int MCMC warmup iterations to send to CmdStan
#' @param refresh  int how often should CmdStan print an update during sampling
#' @param ... additional arguments to pass to CmdStan
#'
#' @return fit a CmdStan fit object
#' @export
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   df <- pn_arx_sim_data(N = 100, ar_X1 = c(0.24), ma_X1 = numeric(0), ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
#'   
#'   df_ppd <- df[91:100,]
#'   fit_pn_arx_model(theta = df$theta[1:90], X = df[1:90,c("X1", "X2")], X_ppd = df_ppd[,c("X1", "X2")])
#' }
fit_pn_arx_model <- function(theta, X, X_ppd, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, ...) {
  stopifnot(is.numeric(theta) && max(abs(theta)) < 2*pi)
  U <- angle_to_unit_vec(theta)
  X_mat <- stats::model.matrix(~., data = as.data.frame(X))
  X_ppd_mat <- stats::model.matrix(~., data = as.data.frame(X_ppd))
  
  model <- instantiate::stan_package_model(
    name = "pn_arx_process",
    package = "pnregstan"
  )
  data_list <- list(N = length(theta),
                    N_ppd = nrow(X_ppd_mat),
                    P = ncol(X_mat),
                    U = U,
                    X = X_mat,
                    X_ppd = X_ppd_mat,
                    # hyperparameters
                    sigma_0 = 100)
  fit <- model$sample(data = data_list, 
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = 2, ...)
  
  return(fit)
}