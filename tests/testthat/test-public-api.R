testthat::test_that("the standardized slab kernels have the expected values", {
  testthat::expect_equal(wendland_kernel(0), 1.5)
  testthat::expect_equal(semicircle_kernel(0), 2 / pi)
  testthat::expect_equal(wendland_kernel(c(-1, 1, 2)), c(0, 0, 0))
  testthat::expect_equal(semicircle_kernel(c(-1, 1, 2)), c(0, 0, 0))
})

testthat::test_that("level_posterior returns normalized posterior probabilities", {
  fit <- level_posterior(
    d = c(-1, 0, 1),
    pi_j = 0.75,
    omega_j = 0.5,
    beta_j = 2,
    sigma = 1,
    likelihood = "gaussian",
    laplace_rate = 0.5,
    quadrature_n = 16L,
    use_exact_wendland = TRUE
  )
  testthat::expect_equal(dim(fit$probability), c(3L, 3L))
  testthat::expect_equal(
    rowSums(fit$probability),
    rep(1, 3),
    tolerance = 1e-10
  )
  testthat::expect_length(fit$posterior_mean, 3L)
})

testthat::test_that("endpoint posterior rules are handled", {
  fit <- level_posterior(
    d = c(-1, 0, 1),
    pi_j = 1,
    omega_j = 0.5,
    beta_j = 2,
    sigma = 1,
    likelihood = "gaussian",
    laplace_rate = 0.5,
    quadrature_n = 16L,
    use_exact_wendland = TRUE
  )
  testthat::expect_equal(fit$posterior_mean, c(0, 0, 0))
  testthat::expect_true(all(fit$probability[, "spike"] == 1))
})
