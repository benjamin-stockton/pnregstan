% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_pgp_geostat_reg_inc.R
\name{create_pgp_geostat_reg_inc_data_list}
\alias{create_pgp_geostat_reg_inc_data_list}
\title{Create Data List for Incomplete Data Projected Gaussian Process Regression}
\usage{
create_pgp_geostat_reg_inc_data_list(loc1, loc2, theta, X)
}
\arguments{
\item{loc1}{Coordinates for the observed sites}

\item{loc2}{Coordinates for the unobserved sites}

\item{theta}{(vector) Incomplete observed angles}

\item{X}{(matrix) Observed predictors; should be complete}
}
\value{
List object to pass to cmdstan
}
\description{
Create Data List for Incomplete Data Projected Gaussian Process Regression
}
\examples{
## Not run:
loc1 <- create_grid_df(10, 10)
loc2 <- create_grid_df(5, 5)
theta <- c(runif(50, 0, 2*pi), rep(NA, 50))
print(theta)
X <- matrix(rnorm(2*100, 0, 1), nrow = 100)

lst <- create_pgp_geostat_reg_inc_data_list(
 loc1 = loc1,
 loc2 = loc2,
 theta = theta,
 X = X)
## End(Not run)
}
