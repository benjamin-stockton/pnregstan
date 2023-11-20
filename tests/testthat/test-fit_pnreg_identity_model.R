test_that("PN Reg Id Sampling Runs", {
  df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
  df_ppd <- pnreg_sim_data(N = 10, mu_X = 0, sigma_X = 1)
  
  fit <- fit_pnreg_identity_model(theta = df$theta,
                                  X = df[,c("X")],
                                  X_ppd = df_ppd[,c("X")],
                                  refresh = 0)
  
  expect_contains(class(fit), "CmdStanFit")
})
