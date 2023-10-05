#' Simulate the effect of a restoration policy over time.
#' @description Using this function, the user can simulate the effect of an
#' intended management policy on the genetic composition of a focal population.
#' The population is assumed to have overlapping generations, and the user can
#' specify two genetic models, either using a simplified average ancestry
#' representation (genetic_model = "point"), or a more detailed model
#' tracking explicit recombination among chromosomes, using genetic_model =
#' "junctions".
#' @inheritParams default_params_doc
#'
#' @param verbose provides verbose output if TRUE.
#' @rawNamespace useDynLib(simRestore)
#' @rawNamespace import(Rcpp)
#' @return tibble with 8 columns: 1) replicate, 2) time (in generations), 3)
#' average frequency of ancestry across all individuals 4) average frequency
#' of ancestry across all males, 5) average frequency of ancestry across all
#' females, 6) number of individuals, 7) number of males and 8) number of
#' females
#' if return_genetics = TRUE, the output is a list containing the above
#' mentioned tibble, called 'results', and a second tibble called 'genetics',
#' with the local ancestry in long format, split out per generation, replicate,
#' individual, sex, linkage group and chromosome (1 or 2). Here, linkage group
#' indicates the focal chromosome (linkage group), and 'chromosome' indicates
#' which of the diploid pair of chromosomes is measured, allowing for phased
#' output if required.
#' @examples
#' sim_pop <- simulate_policy(initial_population_size = 100,
#'                            num_generations = 20,
#'                            starting_freq = 0.2,
#'                            K = 400)
#' plot(sim_pop$results$num_individuals ~ sim_pop$results$t)
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
                            genetic_model = "point",
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
                            ancestry_pull = 1.0,
                            random_mating = FALSE,
                            extra_pair_copulation = 0.0,
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
  if (genetic_model == "junctions" || genetic_model == "junction") {
    use_simplified_model <- FALSE
  }
  if (!(genetic_model %in% c("point", "junctions", "junction"))) {
    stop("Invalid genetic model specified. Did you mean junctions?")
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
                              ancestry_pull,
                              random_mating,
                              extra_pair_copulation,
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
                                       "sex", "linkage_group",
                                       "chromosome", "ancestry")
      } else {
        colnames(output$genetics) <- c("generation", "replicate", "individual",
                                       "sex", "linkage_group", "chromosome",
                                       "position", "ancestry")
      }
      output$genetics <- tibble::as_tibble(output$genetics)
    } else {
      output$genetics <- "extinct"
    }
  }

  return(output)
}
