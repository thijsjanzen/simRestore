
#' Default parameter documentation
#'
#' @param initial_population_size population size at the start
#' @param reproduction_success_rate frequency of females that yield offspring
#' at the end of the breeding season (e.g. a fraction of
#' 1 - reproduction_success_rate of females. This is a joint effect of
#' breeding females getting killed
#' (see \code{female_death_rate}) and other sources of failure to produce
#' offspring. Other sources of failure are calculated from
#' reproduction_success_rate and female_death_rate, such that the resulting
#' reproduction failure rate = 1 - reproduction_success_rate /
#' (1 - female breeding risk);
#' @param reproductive_risk Additional death rate of males and females as a
#' result of breeding (e.g. as a result of protecting the offspring against
#' predators). Provide as a vector where the first index indicates the risk for
#' females, the second the risk for males.
#' @param K carrying capacity
#' @param num_generations number of generations
#' @param pull vector of the number of individuals pulled per year
#' @param put vector of the number of individuals added per year
#' @param starting_freq initial frequency of target ancestry in the population.
#' @param sd_starting_freq variation in initial frequency of target ancestry.
#' @param morgan a vector with the size of each chromosome in morgan, e.g. if
#' a single chromosome is to be simulated a single number will suffice, but
#' for two chromosomes of a size of 1 Morgan, a vector like: c(1, 1) will work.
#' @param establishment_burnin number of generations before establishment
#' @param num_replicates number of replicates
#' @param seed random number seed, if left open, current time is used.
#' @param max_age maximum age an individual can reach.
#' @param mean_number_of_offspring mean number of offspring per female
#' @param sd_number_of_offspring standard deviation of number of offspring per
#' female (assuming the number of offspring is always 0 or larger).
#' @param smin minimum survival rate
#' @param smax maximum survival rate
#' @param b steepness of the survival rate. Negative values indicate a declining
#' survival rate with increasing population size, positive values indicate an
#' increasing survival rate with increasing population size.
#' @param p Density at which the survival rate changes most relative. Expressed
#' in N / K (e.g., for a value of 1.0, the survival rate changes most rapidly
#' around N = K, for a value of 0.5, the survival rate changes most rapidly
#' around N = 0.5K, etc).
#' @param sex_ratio_put sex ratio of individuals that are added (if any) to
#' the population. Sex ratio is expressed as males / (males + females), such
#' that 0.5 indicates an even sex ratio, 0.9 indicates a male biased sex ratio
#' and 0.1 indicates a female biased sex ratio.
#' @param sex_ratio_pull sex ratio of individuals that are removed (if any) from
#' the population. The sex ratio
#' is expressed as males / (males + females), such
#' that 0.5 indicates an even sex ratio, 0.9 indicates a male biased sex ratio
#' and 0.1 indicates a female biased sex ratio.
#' @param sex_ratio_offspring sex ratio of newly born offspring. The sex ratio
#' is expressed as males / (males + females), such
#' that 0.5 indicates an even sex ratio, 0.9 indicates a male biased sex ratio
#' and 0.1 indicates a female biased sex ratio.
#' @param ancestry_put Average ancestry of individuals being used for
#' supplementation. If the target is high focal ancestry (e.g. aiming for
#' focal ancestry of 1.0), ancestry put should reflect this and be set to 1.0 (
#' which is the default value). When supplementing with non-pure individuals,
#' this value can consequently be lowered.
#' @param ancestry_pull Ancestry level below which individuals are allowed to
#' be pulled - this can reduce the effective number of individuals pulled if
#' none of the individuals in the population match this ancestry level. This can
#' be used to selectively only remove those with low target ancestry.
#' @param random_mating by default, simulations assume fixed pair bonding, e.g.
#' each female mates with exactly one male (if available). Alternatively, if
#' random_mating = TRUE, females will mate with a random male, introducing the
#' possibility that some males mate multiple times.
#' @param extra_pair_copulation probability of offspring to be fathered by
#' another male. We assume that all offspring from one mother can have at most
#' two fathers.
#' @param genetic_model The model can either use the point ancestry model
#' ("point") of underlying genetics, which speeds up simulation
#' considerably, but underestimates genetic variation. Alternatively, a more
#' detailed genetic model is available, making use of the theory of junctions,
#' this can be accessed using the option "junctions". Default is "point".
#' @param return_genetics returns a tibble containing all local ancestry
#' information for all individuals. This tibble contains the following
#' informative columns: 1) time (the last generation), 2) replicate,
#' 3) individual, 4) sex (0 = male, 1 = female), 5) Linkage Group
#' (if use_simplified_model == FALSE), 6) chromosome (1 or 2,
#' returning phased results), 7) position (if use_simplified_model == FALSE) and
#' 8) local ancestry (0 or 1).
#' @param verbose provides verbose output if TRUE.
#' @return Nothing
#' @keywords internal
default_params_doc <- function(initial_population_size = 400,
                               reproduction_success_rate = 0.387,
                               breeding_risk = c(0.2, 0.0),
                               K = 400, # nolint
                               num_generations = 20,
                               pull = 0,
                               put = 0,
                               starting_freq = 0.5,
                               sd_starting_freq = 0.05,
                               morgan = c(1),
                               establishment_burnin = 30,
                               num_replicates = 1,
                               seed = NULL,
                               max_age = 6,
                               mean_number_of_offspring = 6,
                               sd_number_of_offspring = 1,
                               smin = 0.5,
                               smax = 0.9,
                               b = -2,
                               p = 0.5,
                               sex_ratio_put = 0.5,
                               sex_ratio_pull = 0.5,
                               sex_ratio_offspring = 0.5,
                               ancestry_put = 1.0,
                               use_simplified_model = TRUE,
                               verbose = TRUE) {
 # Nothing
}
