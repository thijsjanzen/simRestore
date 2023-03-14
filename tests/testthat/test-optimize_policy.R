context("optimize policy")

test_that("simple optimization", {

  vx <- simRestore::optimize_policy(num_generation = 20,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = FALSE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.99)
  testthat::expect_equal(max(vx$results$t), 20)
  testthat::expect_true(is.na(sd(vx$put)))

  vx <- simRestore::optimize_policy(num_generation = 20,
                                    target_frequency = 0.99,
                                    optimize_put = FALSE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)
  testthat::expect_lt(vx$final_freq, 0.99)
  testthat::expect_equal(max(vx$results$t), 20)
  testthat::expect_true(is.na(sd(vx$pull)))

  vx <- simRestore::optimize_policy(num_generation = 20,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.99)
  testthat::expect_equal(max(vx$results$t), 20)
  testthat::expect_true(is.na(sd(vx$pull)))
  testthat::expect_true(is.na(sd(vx$put)))
})
