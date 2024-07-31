test_that("Data list creation is correct", {
  
  loc1 <- create_grid_df(2, 2)
  loc2 <- create_grid_df(5, 5)
  
  dl <- create_pgp_geostat_reg_inc_data_list(
   loc1 = loc1,
   loc2 = loc2,
   theta = c(0, pi/2, NA, 3*pi/2),
   X = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1))
  
  dl_true <- list(N_obs = 3,
                     N_inc = 1,
                     N2 = 25,
                     P = 1,
                     theta1 = c(0, pi/2, 3*pi/2),
                     Xobs = matrix(c(0, 1/4, 3/4), ncol = 1),
                     Xinc = matrix(c(1/2), ncol = 1),
                     loc_obs = loc1[c(1,2,4),],
                     loc_inc = loc1[3,],
                     loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("2nd data list creation is correct", {
  
  loc1 <- create_grid_df(10, 10)
  loc2 <- create_grid_df(5, 5)
  theta <- c(runif(50, 0, 2*pi), rep(NA, 50))
  X <- matrix(rnorm(2*100, 0, 1), nrow = 100)

  dl <- create_pgp_geostat_reg_inc_data_list(
   loc1 = loc1,
   loc2 = loc2,
   theta = theta,
   X = X)
  
  dl_true <- list(N_obs = 50,
                  N_inc = 50,
                  N2 = 25,
                  P = 2,
                  theta1 = theta[1:50],
                  Xobs = X[1:50,],
                  Xinc = X[51:100,],
                  loc_obs = loc1[1:50,],
                  loc_inc = loc1[51:100,],
                  loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("PGP Geostatistical Regression Incomplete Sampling Runs", {
  set.seed(651)

  loc <- create_grid_df(10, 10)
    loc2 <- create_grid_df(5, 5)
    df <- data.frame(X = 1:nrow(loc) / nrow(loc))
    df$theta <- df$X * 2*pi
    ind_miss <- 1:25 * 4
    df$theta_inc <- df$theta
    df$theta_inc[ind_miss] <- NA

    fit <- fit_pgp_geostat_reg_inc_model(loc1 = loc, loc2 = loc2,
                                  theta = df$theta_inc,
                                  X = df$X,
                                  iter_warmup = 10,
                                  iter_sampling = 10,
                                  show_messages = TRUE,
                                  show_exceptions = TRUE)

  expect_contains(class(fit), "CmdStanFit")
})
