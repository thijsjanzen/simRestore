context("simulate simple and nonsimple")

test_that("compare use", {
  vx <- simulate_policy(initial_population_size = 1000,
                        nest_success_rate = 0.387,
                        nesting_risk = 0.2,
                        K = 400,
                        num_generations = 20,
                        pull = 0,
                        put = 0,
                        num_replicates = 100,
                        starting_freq = 0.2,
                        seed = 42,
                        use_simplified_model = FALSE,
                        verbose = FALSE)

  vy <- simulate_policy(initial_population_size = 1000,
                        nest_success_rate = 0.387,
                        nesting_risk = 0.2,
                        K = 400,
                        num_generations = 20,
                        pull = 0,
                        put = 0,
                        num_replicates = 100,
                        starting_freq = 0.2,
                        seed = 42,
                        use_simplified_model = TRUE)

  for (tt in unique(vx$results$t)) {
    if (tt > 1) {
      a <- subset(vx$results, vx$results$t == tt)
      b <- subset(vy$results, vy$results$t == tt)
      vv <- t.test(a$Num_individuals, b$Num_individuals)
      vv2 <- t.test(a$freq_hawaii, b$freq_hawaii)
      testthat::expect_true(vv2$p.value > 0.001)
      testthat::expect_true(vv$p.value > 0.001)
    }
  }
})
