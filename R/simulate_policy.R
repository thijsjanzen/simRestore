#' Simulate the effect of a policy over time.
#' @param initial_population_size population size at the start
#' @param nest_success_rate frequency of nests that yield offspring at the end
#' of the breeding season (e.g. a fraction of 1 - nest_success_rate of nests
#' fail). This is a joint effect of breeding females getting killed
#' (see \code{female_death_rate}) and other sources of failure to complete a
#' nest. Other sources of failure are calculated from nest_success_rate and
#' female_death_rate, such that nest failure rate = 1 - nest_success_rate / (1 -
#' female_death_rate);
#' @param nesting_risk Additional death rate of males and females as a result of
#' protecting the nest. Provide as a vector where the first index indicates
#' the risk for females, the second the risk for males.
#' @param K carrying capacity
#' @param num_generations number of generations
#' @param pull vector of the number of individuals pulled per year
#' @param put vector of the number of individuals added per year
#' @param starting_freq initial hawaii frequency in the population.
#' @param sd_starting_freq variation in initial hawaii frequency.
#' @param morgan size of the chromosome in Morgan
#' @param establishment_burnin number of generations before establishment
#' @param num_replicates number of replicates
#' @param seed random number seed, if left open, current time is used.
#' @param max_age maximum age a duck can reach.
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
#' @param verbose provides verbose output if TRUE.
#' @rawNamespace useDynLib(simRestore)
#' @rawNamespace import(Rcpp)
#' @return tibble
#' @export
simulate_policy <- function(initial_population_size = 400,
                            nest_success_rate = 0.387,
                            nesting_risk = c(0.2, 0.0),
                            K = 400, # nolint
                            num_generations = 20,
                            pull = 0,
                            put = 0,
                            starting_freq = 0.5,
                            sd_starting_freq = 0.05,
                            morgan = 1,
                            establishment_burnin = 30,
                            num_replicates = 1,
                            seed = NULL,
                            max_age = 6,
                            mean_clutch_size = 6,
                            sd_clutch_size = 1,
                            smin = 0.5,
                            smax = 0.9,
                            b = -2,
                            p = 0.5,
                            sex_ratio_put = 0.5,
                            sex_ratio_offspring = 0.5,
                            use_simplified_model = TRUE,
                            verbose = FALSE) {

  shooting <- update_vector(pull, num_generations)
  addition <- update_vector(put, num_generations)

  nest_failure_rate <- 1 - nest_success_rate / (1 - nesting_risk)

  seed <- set_seed(seed)

  if (smin >= smax) {
    stop("smin has to be smaller than smax")
  }

  output <- simulate_complete(initial_population_size,
                              starting_freq,
                              sd_starting_freq,
                              addition,
                              shooting,
                              num_generations,
                              num_replicates,
                              K,
                              morgan,
                              nesting_risk,
                              nest_failure_rate,
                              establishment_burnin,
                              seed,
                              max_age,
                              use_simplified_model,
                              verbose,
                              mean_clutch_size,
                              sd_clutch_size,
                              smin,
                              smax,
                              p,
                              b,
                              sex_ratio_put,
                              sex_ratio_offspring)

  colnames(output$results) <- c("replicate", "t", "freq_hawaii",
                                "freq_hawaii_males", "freq_hawaii_females",
                                "Num_individuals",
                                "Num_males", "Num_females")
  output$results <- tibble::as_tibble(output$results)
  return(output)
}
