#' Optimize putting and/or pulling, where it is assumed that the same amount
#' is applied per generation.
#' @param num_generations number of generations
#' @param target_frequency frequency to aim for
#' @param K carrying capacity
#' @param optimize_pull When set to 0, FALSE or a negative number, it will not
#' be optimized. When negative, the absolute value will be taken as a fixed
#' contribution to each generation (but will not be optimized)
#' @param optimize_put When set to 0, FALSE or a negative number, it will not
#' be optimized. When negative, the absolute value will be taken as a fixed
#' contribution to each generation (but will not be optimized)
#' @param num_replicates Number of replicates per parameter combination to be
#' simulated. Fit of the parameter combination is chosen as the average
#' frequency across replicates. Please note that this does not refer to the
#' number of replicate optimization routine: the number of replicates refers
#' to the number of simulations used PER parameter evaluation in the
#' optimization.
#' @param initial_population_size population size at the start
#' @param nest_success_rate frequency of nests that yield offspring at the end
#' of the breeding season (e.g. a fraction of 1 - nest_success_rate of nests
#' fail). This is a joint effect of breeding females getting killed
#' (see \code{nesting_risk}) and other sources of failure to complete a
#' nest. Other sources of failure are calculated from nest_success_rate and
#' nesting_risk, such that nest failure rate = 1 - nest_success_rate / (1 -
#' nesting_risk);
#' @param nesting_risk Additional death rate of females as a result of
#' protecting the nest.
#' @param num_generations number of generations
#' @param starting_freq initial focal frequency in the population.
#' @param morgan size of the chromosome in Morgan
#' @param establishment_burnin number of generations before establishment
#' @param num_replicates number of replicates
#' @param verbose provide verbose output
#' @param max_age maximum age an individual can reach.
#' @param mean_clutch_size mean number of eggs in a nest
#' @param sd_clutch_size standard deviation of number of eggs in nest (assuming
#' the number of eggs is always 0 or larger).
#' @param smin minimum survival rate
#' @param smax maximum survival rate
#' @param b steepness of the survival rate. Negative values indicate a declining
#' survival rate with increasing population size, positive values indicate an
#' increasing survival rate with increasing population size.
#' @param p Density at which the survival rate changes most relative. Expressed
#' in N / K (e.g., for a value of 1.0, the survival rate changes most rapidly
#' around N = K, for a value of 0.5, the survival rate changes most rapidly
#' around N = 0.5K, etc).
#' @param sex_ratio_put the sex ratio of individuals that are added (if any) to
#' the population. Sex ratio is expressed as males / (males + females), such
#' that 0.5 indicates an even sex ratio, 0.9 indicates a male biased sex ratio
#' and 0.1 indicates a female biased sex ratio.
#' @param sex_ratio_offspring sex ratio of newly born offspring. The sex ratio
#' is expressed as males / (males + females), such
#' that 0.5 indicates an even sex ratio, 0.9 indicates a male biased sex ratio
#' and 0.1 indicates a female biased sex ratio.
#' @param use_simplified_model use a simplified model of underlying genetics?
#' This speeds up simulation considerably, and should be preferred when not
#' interested in high detail genetic changes. For subsequent ADMIXTURE analysis,
#' use_simplified_model should be set to FALSE. Default is TRUE.
#' @return tibble
#' @export
optimize_static <- function(num_generations = 20,
                            target_frequency = 0.99,
                            K = 400, # nolint
                            optimize_put = TRUE,
                            optimize_pull = FALSE,
                            num_replicates = 1,
                            use_simplified_model = TRUE,
                            verbose = FALSE,
                            initial_population_size = 400,
                            nest_success_rate = 0.387,
                            nesting_risk = c(0.2, 0.0),
                            starting_freq = 0.2,
                            morgan = 1,
                            max_age = 6,
                            mean_clutch_size = 6,
                            sd_clutch_size = 1,
                            smin = 0.5,
                            smax = 0.9,
                            b = -2,
                            p = 0.5,
                            sex_ratio_put = 0.5,
                            sex_ratio_offspring = 0.5,
                            establishment_burnin = 30) {

  fit_adding <- function(param,
                         return_results = FALSE) {
    if (param[[1]] < 0) return(Inf)
    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = abs(optimize_pull),
                              put = floor(10^param[[1]]),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_clutch_size = mean_clutch_size,
                              sd_clutch_size = sd_clutch_size,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              use_simplified_model = use_simplified_model,
                              verbose = verbose)

    b <- subset(result$results, result$results$t == num_generations)
    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency) + param[[1]] / 1e6
    if (verbose)  {
      cat(floor(10^param[[1]]), freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return(result$results)
    }
  }

  fit_killing <- function(param,
                          return_results = FALSE) {
    if (param[[1]] < 0) return(100)

    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = param[[1]],
                              put = abs(optimize_put),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_clutch_size = mean_clutch_size,
                              sd_clutch_size = sd_clutch_size,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              use_simplified_model = use_simplified_model,
                              verbose = verbose)


    b <- subset(result$results, result$results$t == num_generations)
    if (length(b$replicate) < 1) {
      return(100)
    }

    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency) + param[[1]] / 1e6
    if (verbose)  {
      cat(param[[1]], freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return(result$results)
    }
  }

  fit_both <- function(param,
                       return_results = FALSE) {
    # normal fit is < 1, this is equal to infinite, but without warning
    if (param[[1]] < 0) return(Inf)
    if (param[[2]] < 0) return(Inf)

    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = (10^param[[1]]),
                              put = (10^param[[2]]),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_clutch_size = mean_clutch_size,
                              sd_clutch_size = sd_clutch_size,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              use_simplified_model = use_simplified_model,
                              verbose = verbose)

    b <- subset(result$results, result$results$t == num_generations)
    if (length(b$replicate) < 1) {
      return(Inf)
    }

    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency) + param[[1]] / 1e6 + param[[2]] / 1e6
    if (verbose)  {
      cat((10^param[[1]]), (10^param[[2]]), freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return(result$results)
    }
  }

  result <- c()

  if (optimize_pull[[1]] > 0 &&
      optimize_put[[1]]  > 0) {

    fit_result <- subplex::subplex(par = c(1, 1), fn = fit_both,
                                   control = list(maxiter = 200,
                                                  reltol = 0.01))
    result$pull <- (10^fit_result$par[[1]])
    result$put  <- (10^fit_result$par[[2]])
    result$results <- fit_both(fit_result$par, return_results = TRUE)

    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     pull = rep(result$pull, num_generations),
                                     put = rep(result$put, num_generations))

    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }
  if (optimize_put > 0 && optimize_pull <= 0) {
    max_val <- log10(K * 2)

    fit_result <- stats::optimize(f = fit_adding,
                                  interval = c(0, max_val))

    while (fit_result$minimum > 0.9 * max_val && fit_result$objective != 1) {
      max_val <- max_val * 2
      fit_result <- stats::optimize(f = fit_adding,
                                    interval = c(0, max_val))
    }

    result$put <- floor(10^fit_result$minimum[[1]])
    result$pull <- abs(optimize_pull)
    result$results <- fit_adding(fit_result$minimum, return_results = TRUE)
    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     put = rep(result$put,
                                               num_generations),
                                     pull = rep(abs(optimize_pull),
                                                num_generations))
    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }
  if (optimize_pull > 0 && optimize_put <= 0) {
    fit_result <-
      stats::optimize(f = fit_killing,
                      interval = c(0, (initial_population_size)),
                      maximum = FALSE)
    start_val <- initial_population_size
    while (fit_result$objective > 1 ||
           fit_result$minimum / start_val > 0.9) {
      # maximum is normally 1, larger indicates failure to converge
      start_val <- start_val * 0.1
      if (start_val < 1) {
        warning("No optimum found, most likely the population is extinct")
        result$pull <- NA
        result$results <- NA
        result$curve <- NA
        result$final_freq <- NA
        return(result)
      }
      fit_result <-
        stats::optimize(f = fit_killing,
                        interval = c(0, start_val),
                        maximum = FALSE)
    }

    result$pull <- floor(fit_result$minimum[[1]])
    result$put  <- abs(optimize_put)
    result$results <- fit_killing(fit_result$minimum, return_results = TRUE)
    if (length(result$results) == 1) {
      warning("No optimum found, most likely the population is extinct")
      result$pull <- NA
      result$results <- NA
      result$curve <- NA
      result$final_freq <- NA
      return(result)
    }

    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     pull = rep(result$pull,
                                                num_generations),
                                     put  = rep(abs(optimize_put),
                                                num_generations))

    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }

  return(result)
}
