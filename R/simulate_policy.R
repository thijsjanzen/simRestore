#' Simulate the effect of a restoration policy over time.
#' @description Using this function, the user can simulate the effect of an
#' intended management policy on the genetic composition of a focal population.
#' The population is assumed to have overlapping generations, and the user can
#' specify two genetic models, either using a simplified average ancestry
#' representation (use_simplified_model = TRUE), or a more detailed model
#' tracking explicit recombination among chromosomes, setting
#' use_simplified_model = FALSE.
#' @inheritParams default_params_doc
#'
#' @param verbose provides verbose output if TRUE.
#' @rawNamespace useDynLib(simRestore)
#' @rawNamespace import(Rcpp)
#' @return tibble
#' @export
simulate_policy <- function(initial_population_size = 400,
                            reproduction_success_rate = 0.387,
                            reproductive_risk = c(0.2, 0.0),
                            K = 400, # nolint
                            num_generations = 20,
                            pull = 0,
                            put = 0,
                            starting_freq = 0.5,
                            sd_starting_freq = 0.05,
                            morgan = c(1),
                            max_age = 6,
                            mean_number_of_offspring = 6,
                            sd_number_of_offspring = 1,
                            genetic_model = "simplified",
                            establishment_burnin = 30,
                            num_replicates = 1,
                            seed = NULL,
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

  shooting <- update_vector(pull, num_generations)
  addition <- update_vector(put, num_generations)

  nest_failure_rate <- 1 - reproduction_success_rate / (1 - reproductive_risk)
  nest_failure_rate <- max(nest_failure_rate)

  seed <- set_seed(seed)

  if (smin >= smax) {
    stop("smin has to be smaller than smax")
  }

  if (length(reproductive_risk) != 2) {
    stop("breeding risk has to be specified for both sexes")
  }

  if (length(morgan) == 1) {
    morgan <- c(morgan)
  }

  use_simplified_model <- TRUE
  if (genetic_model == "junctions") {
    use_simplified_model <- FALSE
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
                              reproductive_risk,
                              nest_failure_rate,
                              establishment_burnin,
                              seed,
                              max_age,
                              use_simplified_model,
                              verbose,
                              mean_number_of_offspring,
                              sd_number_of_offspring,
                              smin,
                              smax,
                              p,
                              b,
                              sex_ratio_put,
                              sex_ratio_pull,
                              sex_ratio_offspring,
                              ancestry_put,
                              return_genetics)



  colnames(output$results) <- c("replicate", "t", "freq_focal_ancestry",
                                "freq_ancestry_males", "freq_ancestry_females",
                                "num_individuals",
                                "num_males", "num_females")
  output$results <- tibble::as_tibble(output$results)

  if (return_genetics) {
    if (ncol(output$genetics) > 0) {
      if (use_simplified_model) {
        colnames(output$genetics) <- c("generation", "replicate", "individual",
                                       "sex", "chromosome", "ancestry")
      } else {
        colnames(output$genetics) <- c("generation", "replicate", "individual",
                                       "sex", "Linkage_Group", "chromosome",
                                       "position", "ancestry")
      }
      output$genetics <- tibble::as_tibble(output$genetics)
    } else {
      output$genetics <- "extinct"
    }
  }

  return(output)
}
