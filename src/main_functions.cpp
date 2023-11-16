// Copyright 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.


#include <thread> // NOLINT [build/c++11]
#include <chrono> // NOLINT [build/c++11]

#include <vector>
#include "./organism.h"
#include "./rand_t.h"
#include "./main_functions.h"
#include "./analysis.h"

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List simulate_complete(int pop_size,
                             float starting_freq,
                             float sd_starting_freq,
                             Rcpp::NumericVector introductions,
                             Rcpp::NumericVector removal,
                             int number_of_generations,
                             int replicates,
                             int K,
                             std::vector<double> morgan,
                             Rcpp::NumericVector nesting_risk,
                             double nest_failure_rate,
                             int establishment_burnin,
                             int seed,
                             int max_age,
                             bool use_simple,
                             bool verbose,
                             double clutch_size_mean,
                             double clutch_size_sd,
                             double smin,
                             double smax,
                             double p,
                             double b,
                             double sex_ratio_put,
                             double sex_ratio_pull,
                             double sex_ratio_offspring,
                             double put_ancestry,
                             double pull_ancestry,
                             bool use_random_mating,
                             double extra_pair_copulation,
                             bool return_genetics) {
  try {
    emp_genome dummy_genome;

    parameters params(pop_size,
                      starting_freq,
                      sd_starting_freq,
                      number_of_generations,
                      K,
                      morgan,
                      nesting_risk[0],
                      nesting_risk[1],
                      nest_failure_rate,
                      establishment_burnin,
                      max_age,
                      dummy_genome,
                      clutch_size_mean,
                      clutch_size_sd,
                      smin,
                      smax,
                      p,
                      b,
                      sex_ratio_put,
                      sex_ratio_pull,
                      sex_ratio_offspring,
                      put_ancestry,
                      pull_ancestry,
                      extra_pair_copulation,
                      use_random_mating);

    std::vector< std::vector< double > > results;

    Rcpp::NumericMatrix base_gen_dummy;
    std::vector< std::vector< double >> anc_inf_dummy;
    Rcpp::NumericVector base_markers;

    Rcpp::NumericMatrix genetic_output;

    if (use_simple) {
      analysis<organism_simple> main_analysis(params,
                                              introductions,
                                              removal,
                                              base_markers,
                                              verbose,
                                              false,
                                              base_gen_dummy,
                                              anc_inf_dummy,
                                              seed);

      for (int r = 0; r < replicates; ++r) {
        output_data run_result = main_analysis.do_analysis();
        for (int t = 0; t < run_result.size(); ++t) {
          std::vector< double > temp =
            {static_cast<double>(run_result[t].replicate + 1),
             static_cast<double>(run_result[t].t + 1),
             run_result[t].frequency,
             run_result[t].frequency_males,
             run_result[t].frequency_females,
             static_cast<double>(run_result[t].pop_size),
             static_cast<double>(run_result[t].num_males),
             static_cast<double>(run_result[t].num_females)};
          results.emplace_back(temp);
        }
        if (return_genetics) {
          main_analysis.export_genetics(&genetic_output);
        }
        main_analysis.increase_replicate();
      }
    } else {
      analysis<organism> main_analysis(params,
                                       introductions,
                                       removal,
                                       base_markers,
                                       verbose,
                                       false,
                                       base_gen_dummy,
                                       anc_inf_dummy,
                                       seed);

      for (int r = 0; r < replicates; ++r) {

        output_data run_result = main_analysis.do_analysis();

        for (int t = 0; t < run_result.size(); ++t) {
          std::vector< double > temp =
            {static_cast<double>(run_result[t].replicate + 1),
             static_cast<double>(run_result[t].t + 1),
             run_result[t].frequency,
             run_result[t].frequency_males,
             run_result[t].frequency_females,
             static_cast<double>(run_result[t].pop_size),
             static_cast<double>(run_result[t].num_males),
             static_cast<double>(run_result[t].num_females)};
          results.push_back(temp);
        }
        if (return_genetics) {
          main_analysis.export_genetics(&genetic_output);
        }

        main_analysis.increase_replicate();
      }
    }

    Rcpp::NumericMatrix output(results.size(), 8);
    for (size_t i = 0; i < results.size(); ++i) {
      output(i, 0) = results[i][0];
      output(i, 1) = results[i][1];
      output(i, 2) = results[i][2];
      output(i, 3) = results[i][3];
      output(i, 4) = results[i][4];
      output(i, 5) = results[i][5];
      output(i, 6) = results[i][6];
      output(i, 7) = results[i][7];
    }

    Rcpp::List output_list;

    if (return_genetics) {
      output_list = Rcpp::List::create(Rcpp::Named("results")  = output,
                         Rcpp::Named("genetics") = genetic_output);
    } else {
      output_list = Rcpp::List::create(Rcpp::Named("results")  = output);
    }
    return(output_list);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}
