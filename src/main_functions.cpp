#include <thread>
#include <chrono>

#include <vector>
#include "Duck.h"
#include "rand_t.h"
#include "main_functions.h"
#include "analysis.h"
#include "genome_tools.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List simulate_complete(int pop_size,
                       float frequency_hawaii_duck,
                       float sd_frequency_hawaii,
                       NumericVector introductions,
                       NumericVector removal,
                       int number_of_generations,
                       int replicates,
                       int K,
                       double morgan,
                       double female_death_rate,
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
                       double sex_ratio_offspring) {

  try {
    if (verbose) {
      Rcout << "performing a simple model: " << use_simple << "\n";
      force_output(); };

    emp_genome dummy_genome;

    parameters params(pop_size,
                      frequency_hawaii_duck,
                      sd_frequency_hawaii,
                      number_of_generations,
                      K,
                      morgan,
                      female_death_rate,
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
                      sex_ratio_offspring);

    std::vector< std::vector< double > > results;

    NumericMatrix base_gen_dummy;
    std::vector< std::vector< double >> anc_inf_dummy;
    NumericVector base_markers;

    if (use_simple) {
      analysis<Duck_simple> main_analysis(params,
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
      }
    } else {
      analysis<Duck> main_analysis(params,
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
      }
    }

    if (verbose) {Rcout << "done simulating, converting to R format\n"; force_output();}
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

    Rcpp::List output_list = List::create( Named("results") = output);
    return(output_list);

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}


void force_output() {
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}

