% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_vm_brms_model.R
\name{fit_vm_brms_model}
\alias{fit_vm_brms_model}
\title{Fit the von Mises Regression Using Stan}
\usage{
fit_vm_brms_model(
  theta,
  X,
  X_ppd,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 500,
  ...
)
}
\arguments{
\item{theta}{vector of angles}

\item{X}{vector or matrix of predictors}

\item{X_ppd}{vector or matrix of holdout predictors}

\item{iter_sampling}{int number of sampling iterations}

\item{iter_warmup}{int number of warmup iterations}

\item{refresh}{int how often to print updates to console}

\item{...}{Other arguments to pass to cmdstan}
}
\value{
CmdStan fit object
}
\description{
Fit the von Mises Regression Using Stan
}
\examples{
if (instantiate::stan_cmdstan_exists()) {
  df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
  fit_vm_brms_model(theta = df$theta, X = df$X, X_ppd = df$X, 
  refresh = 0,
  iter_warmup = 10,
  iter_sampling = 10,
  show_messages = FALSE,
  show_exceptions = FALSE)
}
}
\seealso{
Other models: 
\code{\link{fit_pnreg_gen_model}()},
\code{\link{fit_pnreg_identity_model}()}
}
\concept{models}
