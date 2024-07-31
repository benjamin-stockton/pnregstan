#' Simulate from Projected Gaussian Process Geostatistical Regression Model
#'
#' @param loc Matrix of coordinates of N sites. 
#' @param X Predictor matrix (N x P)
#' @param sigma_w Standard deviation of the first latent angular coordinate
#' @param rho_w Correlation between latent angular coordinates
#' @param rho Length scale parameter of the exponential covariance function
#' @param mu0 Latent bivariate intercept vector
#' @param B Regression coefficient matrix (P x 2)
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return Data frame with simulated PGP, predictors, and location coordinates
#' @export
#'
#' @examples
#' ## Not run:
#' if (FALSE) {
#'   loc <- create_grid_df(10, 10)
#'   loc2 <- create_grid_df(11, 11)
#'   X <- rnorm(100)
#'   B <- matrix(c(1, 2), nrow = 1)
#'  
#'   df <- pgp_geostat_reg_sim_data(loc, X = X, sigma_w = 1, rho_w = 0, rho = 5,
#'          mu0 = c(1,0), B = B)   
#'   
#' }
#' ## End(Not run)
pgp_geostat_reg_sim_data <- function(loc, X, sigma_w, rho_w, rho, mu0, B, iter_sampling = 1000, iter_warmup = 1000, refresh = 0, chains = 2, ...) {
  
  if (!is.matrix(X) & is.numeric(X)) {
    X <- matrix(X, ncol = 1)
  }
  
  stopifnot(nrow(loc) == nrow(X))
  
  data_list <- list(
    N = nrow(loc),
    P = ncol(X),
    loc = loc,
    X = as.matrix(X),
    sigma_w = sigma_w,
    rho_w = rho_w,
    rho = rho,
    mu0 = mu0,
    B = B)
  
  model <- instantiate::stan_package_model(
    name = "pgp_reg_sim",
    package = "pnregstan"
  )
  
  sim <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  theta_sim <- sim$draws(variables = c("theta")) |> posterior::as_draws_matrix()
  
  theta <- theta_sim[sample(iter_sampling, size = 1),] |> as.numeric()
  df <- data.frame(theta = theta %% (2*pi),
                   long = loc$long,
                   lat = loc$lat)
  
  df <- cbind(df, X)
  
  return(df)
}
