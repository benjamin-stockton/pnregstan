#' Fit the Joint Projected Gaussian Process Multivariate Gaussian Process Geostatistical Regression Model to Incomplete Data
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Incomplete observed angles
#' @param Y (vector) (Incomplete) observed responses
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
#'   df$Y <- cos(df$theta) - 0.2 * sin(df$theta) + rnorm(100)
#'   ind_miss <- 1:25 * 4
#'   df$theta_inc[ind_miss] <- NA
#'          
#'   fit_jpgpmgp_geostat_reg_inc_model(
#'   loc1 = loc, loc2 = loc2,
#'   theta = df$theta_inc, 
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
fit_jpgpmgp_geostat_reg_inc_model <- function(loc1, loc2, theta, Y, X, iter_sampling = 1000, iter_warmup = 1000, refresh = 500, chains = 2, ...) {
  
  data_list <- create_jpgpmgp_geostat_reg_inc_data_list(loc1 = loc1, loc2 = loc2, theta = theta, Y = Y, X = X)
  
  model <- instantiate::stan_package_model(
    name = "jpgpmgp_reg_imp2",
    package = "pnregstan"
  )
  
  fit <- model$sample(data = data_list,
                      iter_sampling = iter_sampling,
                      iter_warmup = iter_warmup,
                      refresh = refresh,
                      chains = chains, ...)
  
  return(fit)
}



#' Create Data List for Incomplete Data Joint Projected Gaussian Process Multivariate Gaussian Process Regression
#'
#' @param loc1 Coordinates for the observed sites
#' @param loc2 Coordinates for the unobserved sites
#' @param theta (vector) Incomplete observed angles
#' @param Y (vector) Incomplete responses
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
#' Y <- cos(theta) - 0.2 * sin(theta) + rnorm(100)
#' 
#' lst <- create_jpgpmgp_geostat_reg_inc_data_list(
#'  loc1 = loc1,
#'  loc2 = loc2,
#'  theta = theta,
#'  Y = Y,
#'  X = X)
#' ## End(Not run)
create_jpgpmgp_geostat_reg_inc_data_list <- function(loc1, loc2, theta, Y, X) {
  miss_mat <- matrix(as.numeric(is.na(cbind(theta, Y))),
                     nrow = nrow(loc1))
  n_mis <- apply(miss_mat, 2, sum)
  
  if (n_mis[1] > 0) {
    theta_mis_ind <- which(is.na(theta))
    theta[theta_mis_ind] <- sample(theta[-theta_mis_ind],
                                   size = n_mis[1], replace = TRUE)
  }
  
  if (n_mis[2] > 0) {
    y_mis_ind <- which(is.na(Y))
    Y[y_mis_ind] <- sample(Y[-y_mis_ind],
                           size = n_mis[2], replace = TRUE)
    
  }
  
  if (!is.matrix(X) & is.numeric(X)) {
    X <- matrix(X, ncol = 1)
    P <- 1
  }
  else if (is.matrix(X) && ncol(X) == 1) {
    X <- matrix(X, ncol = 1)
    P <- 1
  } else if (is.matrix(X) && ncol(X) > 1) {
    X <- 
    P <- ncol(X)
  }
  
  stopifnot(is.numeric(theta) && max(abs(theta)) <= 2*pi)
  
  data_list <- list(
    N1 = nrow(loc1),
    N2 = nrow(loc2),
    P = P,
    theta1 = theta,
    Y1 = Y,
    X1 = as.matrix(X),
    R = miss_mat,
    loc1 = loc1,
    loc2 = loc2
  )
  
  return(data_list)
}