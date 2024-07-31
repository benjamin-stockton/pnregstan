#' Simulate from Projected Gaussian Process Geostatistical Model
#'
#' @param N Number of columns of sites
#' @param M Number of rows of sites
#' @param sigma_w Standard deviation of the first latent angular coordinate
#' @param rho_w Correlation between latent angular coordinates
#' @param rho Length scale parameter of the exponential covariance function
#' @param mu Latent bivariate normal mean vector
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return Data frame with simulated PGP and location coordinates
#' @export
#'
#' @examples
#' ## Not run:
#' if (FALSE) {
#'   loc <- create_grid_df(10, 10)
#'   loc2 <- create_grid_df(11, 11)
#'  
#'   df <- pgp_geostat_sim_data(N, M, sigma_w = 1, rho_w = 0, rho = 5, 
#'          mu = c(1,0))   
#'   
#' }
#' ## End(Not run)
pgp_geostat_sim_data <- function(N, M, sigma_w, rho_w, rho, mu, iter_sampling = 1000, iter_warmup = 1000, refresh = 0, chains = 2, ...) {
  loc <- create_grid_df(N, M)
  
  data_list <- list(
    N = N*M,
    loc = loc,
    sigma_w = sigma_w,
    rho_w = rho_w,
    rho = rho,
    mu = mu)
  
  model <- instantiate::stan_package_model(
    name = "pgp_geostat_sim",
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
  
  return(df)
}
