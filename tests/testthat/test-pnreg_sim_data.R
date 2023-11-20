test_that("Dimensions are correct with simulated X", {
  df <- pnreg_sim_data(100, mu_X = 0, sigma_X = 5)
  expect_equal(ncol(df), 4)
})

test_that("Dimensions are correct with passed X", {
  X <- mvtnorm::rmvnorm(100, mean = c(0,0), diag(2))
  B <- matrix(c(1, 0, 
                0, 0.2,
                1, -0.5), 
              byrow = TRUE, ncol = 2)
  
  df <- pnreg_sim_data(100, X = X, B = B)
  expect_equal(ncol(df), 5)
})
