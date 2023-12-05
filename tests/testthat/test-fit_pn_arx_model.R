test_that("PN ARX Sampling Runs", {
  df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
  df_ppd <- pnreg_sim_data(N = 10, mu_X = 0, sigma_X = 1)
  
  fit <- fit_pn_arx_model(theta = df$theta,
                             X = df[,c("X")],
                             X_ppd = df_ppd[,c("X")],
                             refresh = 500,
                          show_messages = TRUE,
                          show_exceptions = FALSE, 
                          chains = 1)
  
  expect_contains(class(fit), "CmdStanFit")
})#288s

test_that("PN ARX Sampling Runs with Fixed X", {
  X <- mvtnorm::rmvnorm(100, mean = c(0,0,0), diag(c(1, 5, 10)))
  X_ppd <- mvtnorm::rmvnorm(10, mean = c(0,0,0), diag(c(1, 5, 10)))
  B <- matrix(c(1, 0, 
                0, 0.2,
                1, -0.5,
                0, 1), 
              byrow = TRUE, ncol = 2)
  
  df <- pnreg_sim_data(N = 100, X = X, B = B)
  
  fit <- fit_pn_arx_model(theta = df$theta,
                             X = X,
                             X_ppd = X_ppd,
                             refresh = 500,
                          show_messages = TRUE,
                          show_exceptions = FALSE)
  
  expect_contains(class(fit), "CmdStanFit")
})
