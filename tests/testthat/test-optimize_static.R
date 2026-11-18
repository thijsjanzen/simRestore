context("optimize policy")

test_that("simple optimization", {
  # Seed and is.finite() check added by Claude (Anthropic AI assistant) (1.1.6): the test was
  # stochastic and failed on CRAN when the pull-only scenario below went
  # extinct, producing final_freq = NaN.
  set.seed(42)
  vx <- simRestore::optimize_static(num_generations = 10,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = FALSE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.95) # lenient
  testthat::expect_equal(max(vx$results$t), 10)
  testthat::expect_true(is.na(sd(vx$put)))

  vx <- simRestore::optimize_static(num_generations = 5,
                                    target_frequency = 0.5,
                                    optimize_put = FALSE,
                                    optimize_pull = TRUE,
                                    initial_population_size = 100,
                                    K = 500,
                                    num_replicates = 1)

  testthat::expect_true(is.finite(vx$final_freq)) # added by Claude (Anthropic AI assistant)
  testthat::expect_lt(vx$final_freq, 0.999)
  testthat::expect_equal(length(vx$results$t), 5)
  testthat::expect_true(is.na(sd(vx$pull)))

  vx <- simRestore::optimize_static(num_generations = 20,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.9)

  testthat::expect_equal(max(vx$results$t), 20)
  testthat::expect_true(is.na(sd(vx$pull)))
  testthat::expect_true(is.na(sd(vx$put)))

  # check verbose output:
  testthat::expect_output(
    vx <- simRestore::optimize_static(num_generations = 3,
                                      target_frequency = 0.5,
                                      optimize_put = TRUE,
                                      optimize_pull = FALSE,
                                      num_replicates = 1,
                                      verbose = TRUE)
  )

  testthat::expect_output(
    vx <- simRestore::optimize_static(num_generations = 3,
                                      target_frequency = 0.5,
                                      optimize_put = FALSE,
                                      optimize_pull = TRUE,
                                      num_replicates = 1,
                                      verbose = TRUE)
  )

  testthat::expect_output(
    vx <- simRestore::optimize_static(num_generations = 3,
                                      target_frequency = 0.7,
                                      optimize_put = TRUE,
                                      optimize_pull = TRUE,
                                      num_replicates = 1,
                                      verbose = TRUE)
  )
})

test_that("fixed optimization", {
  fixed_level <- 10
  vx <- simRestore::optimize_static(num_generations = 10,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = -fixed_level,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.95)
  testthat::expect_equal(max(vx$results$t), 10)
  testthat::expect_equal(vx$pull, fixed_level)

  vx <- simRestore::optimize_static(num_generations = 10,
                                    target_frequency = 0.99,
                                    optimize_put = -fixed_level,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)

  testthat::expect_equal(max(vx$results$t), 10)
  testthat::expect_equal(vx$put, fixed_level)
})

test_that("error", {
  testthat::expect_warning(
    vx <- simRestore::optimize_static(num_generations = 10,
                                      initial_population_size = 3,
                                      reproductive_risk = c(1, 1), # force death
                                      target_frequency = 0.99,
                                      optimize_put = FALSE,
                                      optimize_pull = TRUE,
                                      num_replicates = 1)
  )

  testthat::expect_true(is.na(vx$pull))
  testthat::expect_true(is.na(vx$results))
  testthat::expect_true(is.na(vx$curve))
  testthat::expect_true(is.na(vx$final_freq))
})

test_that("genetics", {
  vx <- simRestore::optimize_static(num_generations = 3,
                                    target_frequency = 0.5,
                                    optimize_put = TRUE,
                                    optimize_pull = FALSE,
                                    num_replicates = 1,
                                    return_genetics = TRUE)
  # just check it is there:
  testthat::expect_true(inherits(vx$genetics, "tbl_df"))

  vx <- simRestore::optimize_static(num_generations = 3,
                                    target_frequency = 0.5,
                                    optimize_put = TRUE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1,
                                    return_genetics = TRUE)
  # just check it is there:
  testthat::expect_true(inherits(vx$genetics, "tbl_df"))
})

# Test written by Claude (Anthropic AI assistant) (1.1.6): the same seed must give identical
# results, and the returned simulation must cover all generations.
test_that("seed reproducibility", {
  run <- function(seed) {
    simRestore::optimize_static(num_generations = 5,
                                target_frequency = 0.5,
                                optimize_put = FALSE,
                                optimize_pull = TRUE,
                                initial_population_size = 100,
                                K = 500,
                                num_replicates = 1,
                                seed = seed)
  }
  vx1 <- run(seed = 7)
  vx2 <- run(seed = 7)
  testthat::expect_equal(vx1$pull, vx2$pull)
  testthat::expect_equal(vx1$final_freq, vx2$final_freq)
  testthat::expect_equal(vx1$results, vx2$results)
  testthat::expect_equal(max(vx1$results$t), 5)
  testthat::expect_true(is.finite(vx1$final_freq))
})
