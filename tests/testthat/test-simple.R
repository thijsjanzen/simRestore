context("simulate point")

test_that("point use", {
  vx <- simulate_policy(initial_population_size = 100,
                        K = 400,
                        num_generations = 10,
                        pull = 0,
                        put = 0,
                        starting_freq = 0.5,
                        seed = 42)

  testthat::expect_equal(max(vx$results$replicate), 1)
  testthat::expect_equal(max(vx$results$t), 10)

  vx <- simulate_policy(initial_population_size = 100,
                        K = 400,
                        num_generations = 20,
                        pull = 10,
                        put = 1000,
                        starting_freq = 0.2,
                        num_replicates = 1,
                        seed = 42)

  testthat::expect_true(tail(vx$results$freq_focal_ancestry, 1) > 0.99)

  # seed NULL test
  # and shooting / addition vector wrong length test
  vx <- simulate_policy(initial_population_size = 100,
                        K = 400,
                        num_generations = 20,
                        pull = c(10, 10),
                        put = c(200, 200),
                        starting_freq = 0.2,
                        num_replicates = 1,
                        seed = NULL)

  testthat::expect_true(tail(vx$results$freq_focal_ancestry, 1) > 0.99)

  # overshooting test
  vx <- simulate_policy(initial_population_size = 200,
                        K = 100,
                        num_generations = 20,
                        pull =  1000,
                        put = 0,
                        starting_freq = 0.2,
                        num_replicates = 1,
                        seed = NULL)

  testthat::expect_true(tail(vx$results$t, 1) < 20)
})
