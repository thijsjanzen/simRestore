#' update input vector
#' @param v vector or single number
#' @param num_generations number of generations
#' @keywords internal
update_vector <- function(v, num_generations) {
  if (length(v) == 1) {
    v <- rep(v, num_generations)
  }

  if (length(v) != num_generations) {
    num_missing <- num_generations - length(v)
    v <- c(v, rep(utils::tail(v, 1), num_missing))
  }

  assertthat::are_equal(length(v), num_generations)

  return(v)
}
#' @keywords internal
set_seed <- function(seed) {
  seed <- sample(x = 1:1e8, size = 1)
  set.seed(seed)
  return(seed)
}
