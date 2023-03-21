get_decay_curve <- function(total_sum, params, num_generations) {
  decay_curve <- stats::dbeta(x = (1:num_generations) / num_generations,
                              shape1 = 10^params[[1]], shape2 = 10^params[[2]])
  decay_curve <- total_sum * decay_curve / sum(decay_curve)
  decay_curve <- round(decay_curve)
  # shouldn't happen:
  decay_curve[decay_curve < 0] <- 0
  decay_curve[is.na(decay_curve)] <- 0
  while (sum(decay_curve) != total_sum) {
    difference <- total_sum - sum(decay_curve)
    index <- sample(seq_along(decay_curve), 1)
    if (difference > 0) decay_curve[index] <- decay_curve[index] + 1
    if (difference < 0) {
      decay_curve[index] <- decay_curve[index] - 1
      decay_curve[index] <- max(decay_curve[index], 0)
    }
    decay_curve[is.na(decay_curve)] <- 0
  }

  return(decay_curve)
}



#' Optimize a policy assuming a fixed total sum across all generations of
#' individuals that can be put or pulled (e.g. a fixed effort). This fixed total
#' sum is distributed across the generations following a beta distribution, and
#' the parameters of this beta distribution are fitted.
#' @param num_generations number of generations
#' @param target_frequency frequency to aim for
#' @param optimize_pull Optimization proceeds such that the sum of all
#' removal is equal to this number. Switch off by setting to zero.
#' @param optimize_put optimization proceeds such that the sum of all
#' addition over all generations is equal to this number. Switch off by setting
#' to zero. The individuals are distributed over time following a beta
#' distribution.
#' @param num_replicates number of replicates per parameter combination to be
#' simulated. Fit of the parameter combination is chosen as the average
#' frequency across replicates.
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
#' @param K carrying capacity
#' @param num_generations number of generations
#' @param starting_freq initial focal ancestry frequency in the population.
#' @param morgan size of the chromosome in Morgan
#' @param establishment_burnin number of generations before establishment
#' @param num_replicates number of replicates
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
#' @param verbose provide verbose output
#' @param use_simplified_model use a simplified model of underlying genetics?
#' This speeds up simulation considerably, and should be preferred when not
#' interested in high detail genetic changes. For subsequent ADMIXTURE analysis,
#' use_simplified_model should be set to FALSE. Default is TRUE.
#' @return tibble
#' @export
optimize_adaptive <- function(num_generations = 20,
                                       target_frequency = 0.99,
                                       optimize_pull = 0,
                                       optimize_put = 0,
                                       num_replicates = 1,
                                       use_simplified_model = TRUE,
                                       verbose = FALSE,
                                       initial_population_size = 400, #K
                                       nest_success_rate = 0.387,
                                       nesting_risk = c(0.2, 0.0),
                                       K = 400, # nolint
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

  optim_adding <- function(param, # function to optimize putting
                           return_results = FALSE) {
    if (min(param) < 0) return(Inf)

    decay_curve <- get_decay_curve(optimize_put, param, num_generations)

    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = optimize_pull,
                              put = decay_curve,
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
                              verbose = FALSE)

    b <- subset(result$results, result$results$t == num_generations)
    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency)
    if (verbose)  {
      cat(10^param[[1]], 10^param[[2]], freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return_val <- c()
      return_val$results <- result$result
      return_val$curve <- decay_curve
      return(return_val)
    }
  }

  optim_adding2 <- function(param,  # function to optimize pulling
                            return_results = FALSE) {
    if (min(param) < 0) return(Inf)
    decay_curve <- get_decay_curve(optimize_pull, param, num_generations)

    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = decay_curve,
                              put = optimize_put,
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
                              verbose = FALSE)

    b <- subset(result$results, result$results$t == num_generations)
    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency)
    if (verbose)  {
      cat(10^param[[1]], 10^param[[2]], freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return_val <- c()
      return_val$results <- result$result
      return_val$curve <- decay_curve
      return(return_val)
    }
  }

  optim_adding3 <- function(param,  # function to optimize pulling and putting
                            return_results = FALSE) {
    if (min(param) < 0) return(Inf)
    decay_curve <- get_decay_curve(optimize_put, c(param[[1]], param[[2]]),
                                   num_generations)

    decay_curve2 <- get_decay_curve(optimize_pull, c(param[[3]], param[[4]]),
                                    num_generations)

    result <- simulate_policy(initial_population_size = initial_population_size,
                              nest_success_rate = nest_success_rate,
                              nesting_risk = nesting_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = decay_curve,
                              put = decay_curve2,
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
                              verbose = FALSE)

    b <- subset(result$results, result$results$t == num_generations)
    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency)
    if (verbose)  {
      cat(10^param[[1]], 10^param[[2]], 10^param[[3]], 10^param[[4]],
          freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return_val <- c()
      return_val$results <- result$result
      return_val$curve <- cbind(decay_curve, decay_curve2)
      return(return_val)
    }
  }

  result <- c()

  if (optimize_put > 0 && optimize_pull == 0) {

    fit_result <- subplex::subplex(par = c(1, 1),
                                   fn = optim_adding,
                                   control = list(reltol = 0.01))

    result$put <- get_decay_curve(optimize_put, fit_result$par, num_generations)
    result$pull <- optimize_pull
    interm_result <- optim_adding(fit_result$par, return_results = TRUE)
    result$results <- interm_result$results
    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     put = interm_result$curve)
    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }

  if (optimize_put == 0 && optimize_pull > 0) {

    fit_result <- subplex::subplex(par = c(1, 1),
                                   fn = optim_adding2,
                                   control = list(reltol = 0.01))

    result$pull <- get_decay_curve(optimize_pull, fit_result$par,
                                   num_generations)
    result$put  <- optimize_put
    interm_result <- optim_adding2(fit_result$par, return_results = TRUE)
    result$results <- interm_result$results
    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     pull = interm_result$curve)
    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }

  if (optimize_put > 0 && optimize_pull > 0) {

    fit_result <- subplex::subplex(par = c(1, 1,
                                           1, 1),
                                   fn = optim_adding3,
                                   control = list(reltol = 0.01))

    result$pull <- get_decay_curve(optimize_pull,
                                   fit_result$par[1:2],
                                   num_generations)
    result$put  <- get_decay_curve(optimize_put,
                                   fit_result$par[3:4],
                                   num_generations)
    interm_result <- optim_adding3(fit_result$par, return_results = TRUE)
    result$results <- interm_result$results
    result$curve   <- tibble::tibble(t = 1:num_generations,
                                     put = interm_result$curve[, 1],
                                     pull = interm_result$curve[, 2])
    result$final_freq <- mean(
      subset(result$results,
             result$results$t == num_generations)$freq_focal_ancestry)
  }

  return(result)
}
