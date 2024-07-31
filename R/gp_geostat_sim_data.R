#' Simulate from a GP Geostatistical Model
#'
#' @param N Number of columns of sites
#' @param M Number of rows of sites
#' @param sigma Standard deviation of the process
#' @param rho Length scale parameter for the exponential covariance function (inverse of spatial decay)
#' @param alpha Magnitude parameter for the exponential covariance function
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return Data frame with simulated GP and location coordinates
#' @export
#'
#' @examples
#' ## Not run:
#' if (instantiate::stan_cmdstan_exists()) {
#'   loc <- create_grid_df(10, 10)
#'   loc2 <- create_grid_df(11, 11)
#'  
#'   df <- gp_geostat_sim_data(N = 10, M = 10, sigma = 1, rho = 10, alpha = 5)   
#'   
#' }
#' ## End(Not run)
gp_geostat_sim_data <- function(N, M, sigma, rho, alpha, iter_sampling = 1000, iter_warmup = 1000, refresh = 0, chains = 2, ...) {
  loc <- create_grid_df(N, M)
  
  data_list <- list(
    N = N*M,
    loc = loc,
    sigma = sigma,
    rho = rho,
    alpha = alpha)
  
  model <- instantiate::stan_package_model(
    name = "geostat_sim_gp",
    package = "pnregstan"
  )
  
  sim <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  z_sim <- sim$draws(variables = c("y")) |> posterior::as_draws_matrix()
  
  z <- z_sim[sample(1:nrow(z_sim), size = 1),] |> as.numeric()
  df <- data.frame(X = z,
                   long = loc$long,
                   lat = loc$lat)
  
  return(df)
}

#' Create Regularly-spaced Grid of Locations
#'
#' @param N Number of columns of sites
#' @param M Number of rows of sites
#'
#' @return data frame with cartesian coordinates (long, lat)
#' @export
#'
#' @examples
#' create_grid_df(10, 10)
create_grid_df <- function(N, M) {
  x <- 1:N
  y <- 1:M
  loc <- expand.grid(long = x, lat = y)[,1:2]
  return(loc)
}