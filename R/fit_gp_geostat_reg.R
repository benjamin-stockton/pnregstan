#' Fit Geostatistical Regression Model with Gaussian Process
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param Y (vector) Observed response
#' @param X (matrix) Observed predictors
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
#'   df <- gp_geostat_sim_data(loc, X, sigma = 1, rho = 2, alpha = 0.5)
#'   
#'   df$Y <- df$X + X
#'   df$X <- X
#'          
#'   fit_gp_geostat_reg_model(
#'   loc1 = loc, loc2 = loc2,
#'   Y = df$Y, 
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
fit_gp_geostat_reg_model <- function(loc1, loc2, Y, X, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_gp_geostat_reg_data_list(loc1 = loc1, loc2 = loc2, Y = Y, X = X)
  
  model <- instantiate::stan_package_model(
    name = "geostat_reg_gp",
    package = "pnregstan"
  )
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  return(fit)
}

#' Create Data List for Complete Data Gaussian Process Regression
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param Y (vector) Observed response
#' @param X (matrix) Observed predictors
#'
#' @return List object to pass to cmdstan
#' @export
#'
#' @examples
#' ## Not run:
#' lst <- create_gp_geostat_reg_data_list(
#'  loc1 = create_grid_df(10, 10),
#'  loc2 = create_grid_df(15, 15),
#'  Y = runif(100, 0, 5),
#'  X = matrix(rnorm(2*100, 0, 1), nrow = 100))
#' ## End(Not run)
create_gp_geostat_reg_data_list <- function(loc1, loc2, Y, X) {
  
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
  
  
  data_list <- list(N1 = nrow(loc1),
                    N2 = nrow(loc2),
                    P = P,
                    Y1 = Y,
                    X1 = X1,
                    loc1 = loc1,
                    loc2 = loc2)
  
  return(data_list)
}