#include <thread>
#include <chrono>

#include <vector>
#include "organism.h"
#include "rand_t.h"
#include "main_functions.h"
#include "analysis.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List simulate_complete(int pop_size,
                       float starting_freq,
                       float sd_starting_freq,
                       NumericVector introductions,
                       NumericVector removal,
                       int number_of_generations,
                       int replicates,
                       int K,
                       std::vector<double> morgan,
                       NumericVector nesting_risk,
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
                       double put_ancestry) {

  try {
    if (verbose) {
      Rcout << "performing a simple model: " << use_simple << "\n";
      force_output(); };

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
                      put_ancestry);

    std::vector< std::vector< double > > results;

    NumericMatrix base_gen_dummy;
    std::vector< std::vector< double >> anc_inf_dummy;
    NumericVector base_markers;

    NumericVector age_at_death;

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
        if (verbose) { Rcpp::Rcout << "starting replicate: " << r << "\n"; }
        output_data run_result = main_analysis.do_analysis();
        for(int t = 0; t < run_result.size(); ++t) {
          std::vector< double > temp = {static_cast<double>(run_result[t].replicate + 1),
                                        static_cast<double>(run_result[t].t + 1),
                                        run_result[t].frequency,
                                        run_result[t].frequency_males,
                                        run_result[t].frequency_females,
                                        static_cast<double>(run_result[t].pop_size),
                                        static_cast<double>(run_result[t].num_males),
                                        static_cast<double>(run_result[t].num_females)};
          results.emplace_back(temp);
        }
        main_analysis.increase_replicate();
        age_at_death = NumericVector(main_analysis.recorded_death_ages.begin(),
                                     main_analysis.recorded_death_ages.end());
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
        if (verbose) { Rcpp::Rcout << "starting replicate: " << r << "\n"; }

        output_data run_result = main_analysis.do_analysis();
        if (verbose) { Rcpp::Rcout << "recording replicate: " << r << "\n"; }

        for(int t = 0; t < run_result.size(); ++t) {

          std::vector< double > temp = {static_cast<double>(run_result[t].replicate + 1),
                                        static_cast<double>(run_result[t].t + 1),
                                        run_result[t].frequency,
                                        run_result[t].frequency_males,
                                        run_result[t].frequency_females,
                                        static_cast<double>(run_result[t].pop_size),
                                        static_cast<double>(run_result[t].num_males),
                                        static_cast<double>(run_result[t].num_females)};
          results.push_back(temp);
        }

        main_analysis.increase_replicate();
        age_at_death = NumericVector(main_analysis.recorded_death_ages.begin(),
                                     main_analysis.recorded_death_ages.end());
      }
    }

    NumericMatrix output(results.size(), 8);
    for(int i = 0; i < results.size(); ++i) {
      output(i, 0) = results[i][0];
      output(i, 1) = results[i][1];
      output(i, 2) = results[i][2];
      output(i, 3) = results[i][3];
      output(i, 4) = results[i][4];
      output(i, 5) = results[i][5];
      output(i, 6) = results[i][6];
      output(i, 7) = results[i][7];
    }

    Rcpp::List output_list = List::create( Named("results") = output,
                                           Named("age_at_death") = age_at_death);
    return(output_list);

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}


void force_output() {
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}

