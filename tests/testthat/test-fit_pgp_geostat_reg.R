test_that("PGP Reg Data list creation is correct", {
  
  loc1 <- create_grid_df(2, 2)
  loc2 <- create_grid_df(5, 5)
  
  dl <- create_pgp_geostat_reg_data_list(
    loc1 = loc1,
    loc2 = loc2,
    theta = c(0, pi/2, pi, 3*pi/2),
    X = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1))
  
  dl_true <- list(N1 = 4,
                  N2 = 25,
                  P = 1,
                  theta1 = c(0, pi/2, pi, 3*pi/2),
                  X1 = matrix(c(0, 1/4, 1/2, 3/4), ncol = 1),
                  loc1 = loc1,
                  loc2 = loc2)
  
  expect_equal(dl, dl_true)
})

test_that("PGP Geostatistical Regression Sampling Runs", {
  set.seed(651)

  loc <- create_grid_df(10, 10)
  loc2 <- create_grid_df(5, 5)
  df <- data.frame(X = 1:nrow(loc) / nrow(loc))
  df$theta <- df$X * 2*pi

  fit <- fit_pgp_geostat_reg_model(loc1 = loc, loc2 = loc2,
                                theta = df$theta,
                                X = df$X,
                                iter_warmup = 10,
                                iter_sampling = 10,
                                show_messages = TRUE,
                                show_exceptions = TRUE)

  expect_contains(class(fit), "CmdStanFit")
})
