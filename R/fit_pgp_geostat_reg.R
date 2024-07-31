#' Fit the Projected Gaussian Process Geostatistical Regression with Complete Data
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta vector of angles
#' @param X matrix of predictors
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
#'   X <- gp_geostat_sim_data(10, 10, 1, 5, 1)$X
#'  
#'   df <- pgp_geostat_reg_sim_data(loc, X, sigma_w = 1, rho_w = 0, rho = 5, 
#'          mu0 = c(1,0), B = matrix(c(-0.5, 1), nrow = 1))
#'          
#'   fit_pgp_geostat_reg_model(
#'   loc1 = loc, loc2 = loc2,
#'   theta = df$theta, 
#'   X = df$X,
#'   refresh = 0,
#'   iter_warmup = 10,
#'   iter_sampling = 10,
#'   show_messages = FALSE,
#'   show_exceptions = FALSE)
#'   
#'   
#' }
#' ## End(Not run)
fit_pgp_geostat_reg_model <- function(loc1, loc2, theta, X, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_pgp_geostat_reg_data_list(loc1 = loc1, loc2 = loc2, theta = theta, X = X)
  
  model <- instantiate::stan_package_model(
    name = "pgp_geostat_reg",
    package = "pnregstan"
  )
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  return(fit)
}



#' Create Data List for Complete Data Projected Gaussian Process Regression
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Observed angles
#' @param X (matrix) Observed predictors
#'
#' @return List object to pass to cmdstan
#' @export
#'
#' @examples
#' ## Not run:
#' lst <- create_pgp_geostat_reg_data_list(
#'  loc1 = create_grid_df(10, 10),
#'  loc2 = create_grid_df(15, 15),
#'  theta = runif(100, 0, 2*pi),
#'  X = matrix(rnorm(2*100, 0, 1), nrow = 100))
#' ## End(Not run)
create_pgp_geostat_reg_data_list <- function(loc1, loc2, theta, X) {
  
  if (!is.matrix(X) & is.numeric(X)) {
    X1 <- matrix(X, ncol = 1)
    P <- 1
  }
  else if (is.matrix(X) && ncol(X) == 1) {
    X1 <- as.matrix(X, ncol = 1)
    P <- 1
  } 
  else if (is.matrix(X)) {
    X1 <- X
    P <- ncol(X)
  }
  
  stopifnot(is.numeric(theta) && max(abs(theta)) <= 2*pi)
  
  data_list <- list(N1 = nrow(loc1),
                    N2 = nrow(loc2),
                    P = P,
                    theta1 = theta,
                    X1 = X1,
                    loc1 = loc1,
                    loc2 = loc2)
  
  return(data_list)
}
