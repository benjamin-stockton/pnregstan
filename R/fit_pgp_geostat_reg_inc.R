#' Fit the Projected Gaussian Process Geostatistical Regression Model to Incomplete Data
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Incomplete observed angles
#' @param X (matrix) Observed predictors; should be complete
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
#'   ind_miss <- 1:25 * 4
#'   df$theta_inc[ind_miss] <- NA
#'          
#'   fit_pgp_geostat_reg_inc_model(
#'   loc1 = loc, loc2 = loc2,
#'   theta = df$theta_inc, 
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
fit_pgp_geostat_reg_inc_model <- function(loc1, loc2, theta, X, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_pgp_geostat_reg_inc_data_list(loc1 = loc1, loc2 = loc2, theta = theta, X = X)
  
  model <- instantiate::stan_package_model(
    name = "pgp_geostat_reg_inc",
    package = "pnregstan"
  )
  print(data_list)
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)

  return(fit)
}



#' Create Data List for Incomplete Data Projected Gaussian Process Regression
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Incomplete observed angles
#' @param X (matrix) Observed predictors; should be complete
#'
#' @return List object to pass to cmdstan
#' @export
#'
#' @examples
#' ## Not run:
#' loc1 <- create_grid_df(10, 10)
#' loc2 <- create_grid_df(5, 5)
#' theta <- c(runif(50, 0, 2*pi), rep(NA, 50))
#' print(theta)
#' X <- matrix(rnorm(2*100, 0, 1), nrow = 100)
#' 
#' lst <- create_pgp_geostat_reg_inc_data_list(
#'  loc1 = loc1,
#'  loc2 = loc2,
#'  theta = theta,
#'  X = X)
#' ## End(Not run)
create_pgp_geostat_reg_inc_data_list <- function(loc1, loc2, theta, X) {
  mis_ind <- which(is.na(theta))
  if (length(mis_ind) > 0) {
    theta_obs <- theta[-mis_ind]
    if (!is.matrix(X) & is.numeric(X)) {
      X_obs <- matrix(X[-mis_ind], ncol = 1)
      X_inc <- matrix(X[mis_ind], ncol = 1)
      P <- 1
    }
    else if (is.matrix(X) && ncol(X) == 1) {
      X_obs <- matrix(X[-mis_ind,], ncol = 1)
      X_inc <- matrix(X[mis_ind,], ncol = 1)
      P <- 1
    } else if (is.matrix(X) && ncol(X) > 1) {
      X_obs <- X[-mis_ind,]
      X_inc <- X[mis_ind,]
      P <- ncol(X)
    }
    loc_obs <- loc1[-mis_ind,]
    loc_inc <- loc1[mis_ind,]
    
  }
  else {
    theta_obs <- theta
    X_obs <- X
    
    stopifnot(length(mis_ind) > 0)
  }
  print(str(theta_obs))
  print(max(abs(theta_obs)))
  stopifnot(is.numeric(theta_obs) && max(abs(theta_obs)) <= 2*pi)
  
  data_list <- list(N_obs = nrow(loc_obs),
                    N_inc = nrow(loc_inc),
                    N2 = nrow(loc2),
                    P = P,
                    theta1 = theta_obs,
                    Xobs = as.matrix(X_obs),
                    Xinc = as.matrix(X_inc), 
                    loc_obs = loc_obs,
                    loc_inc = loc_inc,
                    loc2 = loc2)
  
  return(data_list)
}