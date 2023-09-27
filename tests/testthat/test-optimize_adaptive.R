context("optimize policy decay")

test_that("point optimization", {
  total_put <- 100
  vx <- simRestore::optimize_adaptive(num_generations = 20,
                                              target_frequency = 0.99,
                                              optimize_put = total_put,
                                              optimize_pull = 0,
                                              num_replicates = 1)
 testthat::expect_equal(sum(vx$put), total_put)
 testthat::expect_equal(sum(vx$pull), 0)



 total_pull <- 100
 vx <- simRestore::optimize_adaptive(num_generations = 20,
                                              target_frequency = 0.99,
                                              optimize_put = 0,
                                              optimize_pull = total_pull,
                                              num_replicates = 1)
 testthat::expect_equal(sum(vx$pull), total_pull)
 testthat::expect_equal(sum(vx$put), 0)

 vx <- simRestore::optimize_adaptive(num_generations = 20,
                                              target_frequency = 0.99,
                                              optimize_put = total_put,
                                              optimize_pull = total_pull,
                                              num_replicates = 1)
 testthat::expect_equal(sum(vx$pull), total_pull)
 testthat::expect_equal(sum(vx$put), total_put)


 # test with verbose output:
testthat::expect_output(
  vx <- simRestore::optimize_adaptive(num_generations = 3,
                                              target_frequency = 0.99,
                                              optimize_put = total_put,
                                              optimize_pull = 0,
                                              num_replicates = 1,
                                              verbose = TRUE)
)


testthat::expect_output(
 vx <- simRestore::optimize_adaptive(num_generations = 3,
                                              target_frequency = 0.99,
                                              optimize_put = 0,
                                              optimize_pull = total_pull,
                                              num_replicates = 1,
                                              verbose = TRUE)
)

testthat::expect_output(
 vx <- simRestore::optimize_adaptive(num_generations = 3,
                                              target_frequency = 0.99,
                                              optimize_put = total_put,
                                              optimize_pull = total_pull,
                                              num_replicates = 1,
                                              verbose = TRUE)
)
})

test_that("genetics", {
  vx <- simRestore::optimize_adaptive(num_generations = 5,
                                      target_frequency = 0.99,
                                      optimize_put = 100,
                                      optimize_pull = 0,
                                      num_replicates = 1,
                                      return_genetics = TRUE)
  # just check it is there:
  testthat::expect_true(inherits(vx$genetics, "tbl_df"))

  vx <- simRestore::optimize_adaptive(num_generations = 5,
                                      target_frequency = 0.99,
                                      optimize_put = 0,
                                      optimize_pull = 100,
                                      num_replicates = 1,
                                      return_genetics = TRUE)
  # just check it is there:
  testthat::expect_true(inherits(vx$genetics, "tbl_df"))

  vx <- simRestore::optimize_adaptive(num_generations = 5,
                                      target_frequency = 0.99,
                                      optimize_put = 100,
                                      optimize_pull = 30,
                                      num_replicates = 1,
                                      return_genetics = TRUE)
  # just check it is there:
  testthat::expect_true(inherits(vx$genetics, "tbl_df"))
})
