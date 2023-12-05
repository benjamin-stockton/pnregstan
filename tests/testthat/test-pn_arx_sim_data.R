test_that("Correct dimensions", {
  df <- pn_arx_sim_data(N = 100, ar_X1 = c(0.24), ma_X1 = numeric(0), ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
  expect_equal(nrow(df), 100)
  expect_equal(ncol(df), 5)
})

test_that("Verify colnames", {
  df <- pn_arx_sim_data(N = 100, ar_X1 = c(0.24), ma_X1 = numeric(0), ar_X2 = c(0.75), ma_X2 = c(-0.25, 0.25))
  cnames <- colnames(df)
  expect_equal(cnames, c("theta", "U1", "U2", "X1", "X2"))
})
