#' Fit the Projected Normal ARX Process with Incomplete Angular Data
#'
#' @param theta vector of angles
#' @param X Matrix of predictors
#' @param X_ppd Matrix of hold-out predictors
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return CmdStan fit
#' @export
#'
#' @examples
#' ## Not run:
#' if (instantiate::stan_cmdstan_exists()) {
#'   df <- pn_arx_sim_data(N = 100, 
#'   ar_X1 = c(0.24), ma_X1 = numeric(0), 
#'   ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
#'   s <- sample(1:90, size = 10)
#'   df[s, "theta"] <- NA
#'   
#'   df_ppd <- df[91:100,]
#'   fit_pn_arx_id_inc_model(theta = df$theta[1:90], 
#'   X = df[1:90,c("X1", "X2")], 
#'   X_ppd = df_ppd[,c("X1", "X2")],
#'   refresh = 0,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   show_messages = FALSE,
#'   show_exceptions = FALSE)
#'   
#'   
#' }
#' ## End(Not run)
fit_pn_arx_id_inc_model <- function(theta, X, X_ppd, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_pnarxid_inc_data_list(theta = theta, X = X, X_ppd = X_ppd, sigma_0 = 100)
  
  model <- instantiate::stan_package_model(
    name = "pn_arx_id_inc",
    package = "pnregstan"
  )
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)

  return(fit)
}

#' Create Data List for PN ARXID Incomplete Model
#'
#' @param theta vector of angles
#' @param X Matrix of predictors
#' @param X_ppd Matrix of hold-out predictors
#' @param sigma_0 Hyperparameter for mu_0
#'
#' @return data_list a list to pass to Stan
#' @export
#'
#' @examples
#' df <- pn_arx_sim_data(N = 100, 
#'   ar_X1 = c(0.24), ma_X1 = numeric(0), 
#'   ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
#'   s <- sample(1:90, size = 10)
#'   df[s, "theta"] <- NA
#' create_pnarxid_inc_data_list(df$theta[1:90], 
#'     df[1:90, c("X1", "X2")], 
#'     X_ppd = df[91:100, c("X1", "X2")])
create_pnarxid_inc_data_list <- function(theta, X, X_ppd, sigma_0 = 100) {
  theta_obs <- theta[which(!is.na(theta))]
  stopifnot(is.numeric(theta_obs) && max(abs(theta_obs)) <= 2*pi)
  U <- angle_to_unit_vec(theta_obs)
  data_list <- list(N_obs = sum(!is.na(theta)),
                    N_mis = sum(is.na(theta)),
                    mis_ind = which(is.na(theta)),
                    obs_ind = which(!is.na(theta)),
                    N_ppd = nrow(X_ppd),
                    K = ncol(X),
                    U_obs = U,
                    X = X,
                    X_ppd = X_ppd, 
                    sigma_0 = 100)
  
  return(data_list)
}
