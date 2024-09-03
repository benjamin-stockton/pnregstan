test_that("Data list creation is correct", {
  set.seed(123)
  loc1 <- create_grid_df(2, 2)
  loc2 <- create_grid_df(5, 5)
  
  dl <- create_jpgpmgp_geostat_reg_inc_data_list(
    loc1 = loc1,
    loc2 = loc2,
    theta = c(0, pi/2, pi, 3*pi/2),
    Y = 1:4,
    X = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1))
  
  dl_true <- list(N1 = 4,
                  N2 = 25,
                  P = 1,
                  theta1 = c(0, pi/2, pi, 3*pi/2),
                  Y1 = 1:4,
                  X1 = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1),
                  R = matrix(rep(0, 8), ncol = 2),
                  loc1 = loc1,
                  loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("2nd data list creation is correct", {
  set.seed(123)
  loc1 <- create_grid_df(2, 2)
  loc2 <- create_grid_df(5, 5)
  
  dl <- create_jpgpmgp_geostat_reg_inc_data_list(
    loc1 = loc1,
    loc2 = loc2,
    theta = c(0, pi/2, NA, 3*pi/2),
    Y = 1:4,
    X = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1))
  
  dl_true <- list(N1 = 4,
                  N2 = 25,
                  P = 1,
                  theta1 = c(0, pi/2, 3*pi / 2, 3*pi/2),
                  Y1 = 1:4,
                  X1 = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1),
                  R = matrix(c(0,0,1,rep(0, 5)), ncol = 2),
                  loc1 = loc1,
                  loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("3rd data list creation is correct", {
  set.seed(123)
  loc1 <- create_grid_df(2, 2)
  loc2 <- create_grid_df(5, 5)
  
  X <- matrix(c(0, 1/4, 1/2, 3/4,
                1, 3, 4, 1), ncol = 2)
  
  dl <- create_jpgpmgp_geostat_reg_inc_data_list(
    loc1 = loc1,
    loc2 = loc2,
    theta = c(0, pi/2, NA, 3*pi/2),
    Y = 1:4,
    X = X)
  
  dl_true <- list(N1 = 4,
                  N2 = 25,
                  P = 2,
                  theta1 = c(0, pi/2, 3*pi / 2, 3*pi/2),
                  Y1 = 1:4,
                  X1 = X, 
                  R = matrix(c(0,0,1,rep(0, 5)), ncol = 2),
                  loc1 = loc1,
                  loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("JPGPMGP Geostatistical Regression Incomplete Sampling Runs", {
  set.seed(651)
  
  loc <- create_grid_df(6, 6)
  loc2 <- create_grid_df(5, 5)
  df <- data.frame(X = 1:nrow(loc) / nrow(loc))
  df$theta <- df$X * 2*pi
  df$Y <- cos(df$theta) - 0.2 * sin(df$theta) + rnorm(nrow(loc))
  ind_miss <- 1:6 * 6
  df$theta_inc <- df$theta
  df$theta_inc[ind_miss] <- NA
  
  fit <- fit_jpgpmgp_geostat_reg_inc_model(loc1 = loc, loc2 = loc2,
                                       theta = df$theta_inc,
                                       Y = df$Y,
                                       X = df$X,
                                       iter_warmup = 5,
                                       iter_sampling = 5,
                                       show_messages = FALSE,
                                       show_exceptions = FALSE)
  
  expect_contains(class(fit), "CmdStanFit")
})
