test_that("vM Reg Sampling Runs", {
  df <- pnreg_sim_data(N = 100, mu_X = 0, sigma_X = 1)
  df_ppd <- pnreg_sim_data(N = 10, mu_X = 0, sigma_X = 1)
  
  fit <- fit_vm_brms_model(theta = df$theta,
                                  X = df[,c("X")],
                                  X_ppd = df_ppd[,c("X")],
                                  refresh = 0,
                           iter_warmup = 1000,
                           iter_sampling = 1000,
                           show_exceptions = FALSE,
                           show_messages = FALSE)
  
  expect_contains(class(fit), "CmdStanFit")
})
