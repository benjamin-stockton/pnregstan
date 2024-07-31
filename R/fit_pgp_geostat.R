#' Fit the Projected Gaussian Process Geostatistical Model with Complete Data
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta vector of angles
#' @param iter_sampling int iterations for sampling
#' @param iter_warmup int iterations for warmup
#' @param refresh int how often to print to screen
#' @param chains int number of MCMC chains
#' @param ... other arguments to pass to CmdStan
#'
#' @return Cmdstan fit object
#' @export
#'
#' @examples
#' ## Not run:
#' if (FALSE) {
#'   loc <- create_grid_df(10, 10)
#'   loc2 <- create_grid_df(11, 11)
#'  
#'   df <- pgp_geostat_sim_data(loc, sigma_w = 1, rho_w = 0, rho = 5, 
#'          mu = c(1,0))
#'          
#'   fit_pgp_geostat_model(
#'   loc1 = loc, loc2 = loc2,
#'   theta = df$theta,
#'   refresh = 0,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   show_messages = FALSE,
#'   show_exceptions = FALSE)
#'   
#'   
#' }
#' ## End(Not run)
fit_pgp_geostat_model <- function(loc1, loc2, theta, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_pgp_geostat_reg_data_list(loc1 = loc1, loc2 = loc2, theta = theta)
  
  model <- instantiate::stan_package_model(
    name = "wang_gelfand_pgp",
    package = "pnregstan"
  )
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  return(fit)
}



#' Create Data List for Complete Data Projected Gaussian Process Geostatistical Model
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Observed angles
#'
#' @return List object to pass to cmdstan
#' @export
#'
#' @examples
#' ## Not run:
#' lst <- create_pgp_geostat_data_list(
#'  loc1 = create_grid_df(10, 10),
#'  loc2 = create_grid_df(15, 15),
#'  theta = runif(100, 0, 2*pi))
#' ## End(Not run)
create_pgp_geostat_data_list <- function(loc1, loc2, theta) {
  theta_obs <- theta[which(!is.na(theta))]
  loc1 <- loc1[which(!is.na(theta)),]
  stopifnot(is.numeric(theta_obs) && max(abs(theta_obs)) <= 2*pi)
  
  data_list <- list(N1 = nrow(loc1),
                    N2 = nrow(loc2),
                    theta1 = theta_obs,
                    loc1 = loc1,
                    loc2 = loc2)
  
  return(data_list)
}
