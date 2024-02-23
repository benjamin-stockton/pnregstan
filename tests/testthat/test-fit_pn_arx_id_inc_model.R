test_that("Data list creation is correct", {
  set.seed(123)
  df <- pn_arx_sim_data(N = 100,
                        ar_X1 = c(0.24), ma_X1 = numeric(0),
                        ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
  s <- sample(1:90, size = 10)
  df[s, "theta"] <- NA
  
  theta_obs <- df$theta[intersect(which(!is.na(df$theta)), 1:90)]
  U <- angle_to_unit_vec(theta_obs)
  dl <- create_pnarxid_inc_data_list(df$theta[1:90], df[1:90, c("X1", "X2")], X_ppd = df[91:100, c("X1", "X2")])
  
  dl_true <- list(
    N_obs = 80,
    N_mis = 10,
    mis_ind = which(is.na(df$theta)),
    obs_ind = setdiff(1:90, s),
    N_ppd = 10,
    K = 2,
    U_obs = U,
    X = df[1:90,c("X1", "X2")],
    X_ppd = df[91:100,c("X1", "X2")], 
    sigma_0 = 100
  )
  
  expect_equal(dl, dl_true)
})

test_that("PN ARXID Incomplete Sampling Runs", {
  df <- pn_arx_sim_data(N = 100,
                        ar_X1 = c(0.24), ma_X1 = numeric(0),
                        ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
  s <- sample(1:90, size = 10)
  df[s, "theta"] <- NA
  
  df_ppd <- df[91:100,]
  df <- df[1:90,]
  fit <- fit_pn_arx_id_inc_model(theta = df$theta,
                          X = df[,c("X1", "X2")],
                          X_ppd = df_ppd[,c("X1", "X2")],
                          iter_sampling = 10,
                          iter_warmup = 10,
                          refresh = 0,
                          show_messages = TRUE,
                          show_exceptions = FALSE)
  
  expect_contains(class(fit), "CmdStanFit")
})
