test_that("Single angle conversion works", {
  U <- angle_to_unit_vec(pi)
  
  expect_equal(sum(U^2), 1)
})

test_that("Multiple angle conversion works", {
  theta <- seq(0, 2*pi-0.1, length.out = 10)
  U <- angle_to_unit_vec(theta)
  U_sum <- sqrt(U[,1]^2 + U[,2]^2)
  
  expect_equal(U_sum, rep(1, 10))
})
