

testthat::test_that("in case of the 1-dimensional VMD the output is valid!", {

  alpha = 2000
  tau = 0
  K = 4              # different compared to the examples section
  DC = FALSE
  init = 0           # different compared to the examples section
  tol = 1e-6

  set.seed(1)
  res_1d = vmd(data = f_sig,
               alpha = alpha,
               tau = tau,
               K = K,
               DC = DC,
               init = init,
               tol = tol,
               verbose = FALSE)

  testthat::expect_true( inherits(res_1d, 'list') & length(res_1d) == 3 & all(dim(res_1d$u) == c(750, 4)) & all(dim(res_1d$u_hat) == c(750, 4)) & inherits(res_1d$omega, 'matrix') )
})


testthat::test_that("in case of the 2-dimensional VMD the output is valid!", {

  alpha = 5000
  tau = 0.25
  K = 2
  DC = TRUE
  init = 0         # different compared to the examples section
  tol = 1e-7

  set.seed(2)
  res_2d = vmd(data = data,
               alpha = alpha,
               tau = tau,
               K = K,
               DC = DC,
               init = init,
               tol = tol,
               verbose = FALSE)

  testthat::expect_true( inherits(res_2d, 'list') & length(res_2d) == 3 & all(dim(res_2d$u) == c(20, 20, 2)) & all(dim(res_2d$u_hat) == c(20, 20, 2)) & inherits(res_2d$omega, 'array') )
})

