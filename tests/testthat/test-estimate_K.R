
testthat::test_that("the 'estimate_K' function returns the expected output!", {

  data("arrhythmia")

  default_vmd_params = list(alpha = 500,
                            tau = 0,
                            DC = FALSE,
                            init = 1,
                            tol = 1e-6)

  res_k = estimate_k_modes(signal_1d = arrhythmia[['MLII']],
                           cor_thresh = 0.5,                             # increase the correlation coefficient to return the results faster
                           default_vmd_params = default_vmd_params,
                           min_K = 2,
                           seed = 1,
                           verbose = FALSE)

  testthat::expect_true( inherits(res_k, 'numeric') & res_k == 2 )
})
