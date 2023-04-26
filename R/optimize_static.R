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
#'
#' @inheritParams default_params_doc
#'
#' @return tibble
#' @export
optimize_static <- function(target_frequency = 0.99,
                            initial_population_size = 400,
                            reproduction_success_rate = 0.387,
                            reproductive_risk = c(0.2, 0.0),
                            K = 400, # nolint
                            num_generations = 20,
                            optimize_put = TRUE,
                            optimize_pull = FALSE,
                            starting_freq = 0.2,
                            sd_starting_freq = 0.05,
                            morgan = c(1),
                            establishment_burnin = 30,
                            num_replicates = 1,
                            max_age = 6,
                            mean_number_of_offspring = 6,
                            sd_number_of_offspring = 1,
                            genetic_model = "simplified",
                            smin = 0.5,
                            smax = 0.9,
                            b = -2,
                            p = 0.5,
                            sex_ratio_put = 0.5,
                            sex_ratio_pull = 0.5,
                            sex_ratio_offspring = 0.5,
                            ancestry_put = 1.0,
                            verbose = FALSE,
                            return_genetics = FALSE) {


  fit_adding <- function(param,
                         return_results = FALSE,
                         return_gens = FALSE) {
    result <- simulate_policy(initial_population_size = initial_population_size,
                              reproduction_success_rate =
                                reproduction_success_rate,
                              reproductive_risk = reproductive_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = abs(optimize_pull),
                              put = floor(10^param[[1]]),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_number_of_offspring =
                                mean_number_of_offspring,
                              sd_number_of_offspring = sd_number_of_offspring,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              genetic_model = genetic_model,
                              verbose = verbose,
                              return_genetics = return_gens)

    b <- subset(result$results, result$results$t == num_generations)
    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency) + param[[1]] / 1e6
    if (verbose)  {
      cat(floor(10^param[[1]]), freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return(result)
    }
  }

  fit_killing <- function(param,
                          return_results = FALSE,
                          return_gen = FALSE) {

    result <- simulate_policy(initial_population_size = initial_population_size,
                              reproduction_success_rate =
                                reproduction_success_rate,
                              reproductive_risk = reproductive_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = param[[1]],
                              put = abs(optimize_put),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_number_of_offspring =
                                mean_number_of_offspring,
                              sd_number_of_offspring = sd_number_of_offspring,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              genetic_model = genetic_model,
                              verbose = verbose,
                              return_genetics = return_gen)

    b <- subset(result$results, result$results$t == num_generations)
    if (length(b$replicate) < 1) {
      if (!return_results) {
        return(100)
      } else {
        return(result)
      }
    }

    freq <- mean(b$freq_focal_ancestry)
    fit <- abs(freq - target_frequency) + param[[1]] / 1e6
    if (verbose)  {
      cat(param[[1]], freq, fit, "\n")
    }
    if (!return_results) {
      return(fit)
    } else {
      return(result)
    }
  }

  fit_both <- function(param,
                       return_results = FALSE,
                       return_gen = FALSE) {
    # normal fit is < 1, this is equal to infinite, but without warning
    if (param[[1]] < 0) return(Inf)
    if (param[[2]] < 0) return(Inf)

    result <- simulate_policy(initial_population_size =
                                initial_population_size,
                              reproduction_success_rate =
                                reproduction_success_rate,
                              reproductive_risk = reproductive_risk,
                              K = K,
                              num_generations = num_generations,
                              pull = (10^param[[1]]),
                              put = (10^param[[2]]),
                              starting_freq = starting_freq,
                              morgan = morgan,
                              establishment_burnin = establishment_burnin,
                              num_replicates = num_replicates,
                              max_age = max_age,
                              mean_number_of_offspring =
                                mean_number_of_offspring,
                              sd_number_of_offspring =
                                sd_number_of_offspring,
                              smin = smin,
                              smax = smax,
                              b = b,
                              p = p,
                              sex_ratio_put = sex_ratio_put,
                              sex_ratio_offspring = sex_ratio_offspring,
                              genetic_model = genetic_model,
                              verbose = verbose,
                              return_genetics = return_gen)

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
      return(result)
    }
  }




  result <- c()

  if (optimize_pull > 0 &&
      optimize_put  > 0) {

    fit_result <- subplex::subplex(par = c(1, 1), fn = fit_both,
                                   control = list(maxiter = 200,
                                                  reltol = 0.01))
    result$pull <- (10^fit_result$par[[1]])
    result$put  <- (10^fit_result$par[[2]])

    all_res <- fit_both(fit_result$par, return_results = TRUE,
                        return_gen = return_genetics)

    result$results <- all_res$results
    if (return_genetics) result$genetics <- all_res$genetics

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

    all_res <- fit_adding(fit_result$minimum,
                          return_results = TRUE,
                          return_gens = return_genetics)

    result$results <- all_res$results
    if (return_genetics) result$genetics <- all_res$genetics
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
    fit_result <- stats::optimize(f = fit_killing,
                                  interval = c(0, (initial_population_size)),
                                  maximum = FALSE)
    start_val <- initial_population_size
    while (fit_result$objective > 1) {
      # maximum is normally 1, larger indicates failure to converge
      start_val <- start_val * 0.95
      fit_result <- stats::optimize(f = fit_killing,
                                    interval = c(0, start_val),
                                    maximum = FALSE)
      if (start_val < 1) {
        warning("No optimum found, most likely pulling has no effect")
        result$pull <- NA
        result$results <- NA
        result$curve <- NA
        result$final_freq <- NA
        return(result)
      }
    }

    result$pull <- floor(fit_result$minimum[[1]])
    result$put  <- abs(optimize_put)
    all_res <- fit_killing(fit_result$minimum, return_results = TRUE,
                           return_gen = return_genetics)


    if (return_genetics) result$genetics <- all_res$genetics
    if (length(result$results) == 1) {
      warning("No optimum found, most likely the population is extinct")
      result$pull <- NA
      result$results <- NA
      result$curve <- NA
      result$final_freq <- NA
      return(result)
    } else {
      result$results <- all_res$results
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
