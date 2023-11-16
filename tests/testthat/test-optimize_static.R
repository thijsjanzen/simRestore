context("optimize policy")

test_that("simple optimization", {
  vx <- simRestore::optimize_static(num_generations = 10,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = FALSE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.95) # lenient
  testthat::expect_equal(max(vx$results$t), 10)
  testthat::expect_true(is.na(sd(vx$put)))

  vx <- simRestore::optimize_static(num_generations = 3,
                                    target_frequency = 0.5,
                                    optimize_put = FALSE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)

  testthat::expect_lt(vx$final_freq, 0.999)
  testthat::expect_equal(max(vx$results$t), 3)
  testthat::expect_true(is.na(sd(vx$pull)))

  vx <- simRestore::optimize_static(num_generations = 10,
                                    target_frequency = 0.99,
                                    optimize_put = TRUE,
                                    optimize_pull = TRUE,
                                    num_replicates = 1)
  testthat::expect_gt(vx$final_freq, 0.95)
  testthat::expect_equal(max(vx$results$t), 10)
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
                                    target_frequency = 0.2,
                                    optimize_put = FALSE,
                                    optimize_pull = TRUE,
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
