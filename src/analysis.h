#ifndef analysis_hpp
#define analysis_hpp

#include <vector>
#include <algorithm>

#include "rand_t.h"
#include "main_functions.h"
#include "genome_tools.h"

template <typename> struct tag { };

struct parameters {

  parameters() {
    pop_size = 100;
    frequency_hawaii_duck = 0.2;
    sd_frequency_hawaii = 0.05;
    number_of_generations = 20;
    K = 400;
    morgan = 1.0;
    female_death_rate = 0.2;
    nest_failure_rate = 0.387;
    establishment_burnin = 50;
    max_age = 6;
    clutch_size_mean = 6.0;
    clutch_size_sd = 1.0;
    smin = 0.5;
    smax = 0.9;
    p = 1.0;
    b = -2.0;
    sex_ratio_put = 0.5;
    sex_ratio_offspring = 0.5;
  }

  parameters(int p,
             float f,
             float s_f,
             int num_gen,
             int k,
             double m,
             double f_d_r,
             double n_f_r,
             int burnin,
             int m_a,
             emp_genome& emp_gen_in,
             double cs_m,
             double cs_sd,
             double smin_in,
             double smax_in,
             double p_in,
             double b_in,
             double sex_ratio,
             double sex_ratio_o) :
  pop_size(p),
  frequency_hawaii_duck(f),
  sd_frequency_hawaii(s_f),
  number_of_generations(num_gen),
  K(k),
  morgan(m),
  female_death_rate(f_d_r),
  nest_failure_rate(n_f_r),
  establishment_burnin(burnin),
  max_age(m_a),
  clutch_size_mean(cs_m),
  clutch_size_sd(cs_sd),
  smin(smin_in),
  smax(smax_in),
  p(p_in),
  b(b_in),
  sex_ratio_put(sex_ratio),
  sex_ratio_offspring(sex_ratio_o),
  empgen(emp_gen_in){
  }

  int pop_size;
  float frequency_hawaii_duck;
  float sd_frequency_hawaii;
  int number_of_generations;
  int K;
  double morgan;
  double female_death_rate;
  double nest_failure_rate;
  int establishment_burnin;
  int max_age;

  double clutch_size_mean;
  double clutch_size_sd;

  double smin;
  double smax;
  double p;
  double b;

  double sex_ratio_put;
  double sex_ratio_offspring;

  emp_genome empgen; // only used in molecular analysis
};

template <typename BIRD>
class analysis {
public:
  analysis() {
    params = parameters();
  }

  analysis(const parameters& param_in,
           const NumericVector& put,
           const NumericVector& take,
           const NumericVector& mark,
           bool verbose_output,
           bool molecular_data_flag,
           const NumericMatrix& base_gen,
           const std::vector< std::vector< double >>& anc_inf,
           int seed) :
  params(param_in),
  introductions(put),
  removal(take),
  markers(mark),
  verbose(verbose_output),
  using_molecular_data(molecular_data_flag),
  base_genomes(base_gen),
  anc_info(anc_inf) {
    rndgen = rnd_t(seed);
    replicate = 0;
  }

  parameters params;

  int r;

  NumericVector introductions;
  NumericVector removal;
  NumericVector markers;

  bool verbose;
  rnd_t rndgen;
  bool using_molecular_data;
  NumericMatrix base_genomes;
  std::vector< std::vector< double >> anc_info;
  std::vector< Duck_emp > pure;
  std::vector< BIRD > output_pop;

  // member functions
  void increase_replicate() {
    replicate++;
  }

  void set_pure(const std::vector<Duck_emp>& p) {
    pure = p;
  }

  Duck_emp draw_pure() {
    int index = rndgen.random_number(static_cast<int>(pure.size()));
    return pure[index];
  }

  output_data do_analysis() {
    std::vector< BIRD > base_pop = create_base_pop(tag<BIRD>{});

    output_data run_result = simulate_policy(base_pop);

    return run_result;
  }

private:
  int replicate;

  output_data simulate_policy(const std::vector< BIRD >& base_pop) {
    std::vector< BIRD > males;
    std::vector< BIRD > females;

    int num_females = static_cast<int>((base_pop.size() / 2.0) * 1.0);

    for(int i = 0; i < num_females; ++i) {
      females.push_back(base_pop[i]);
    }
    for(int i = num_females; i < base_pop.size(); ++i) {
      males.push_back(base_pop[i]);
    }

    output_data frequencies;
    std::array<double, 3> f1 = calc_freq_hawaii<BIRD>(females, males);

    frequencies.add_slice(replicate, 0, f1,
                          static_cast<int>(females.size() + males.size()),
                          males.size(),
                          females.size());

    for (int t = 1; t < params.number_of_generations; ++t) {

      update_pop(females,
                 males,
                 introductions[t],
                 removal[t]);

      if(males.empty() && females.empty()) {
        if (verbose) {
           Rcpp::Rcout << "population went extinct\n"; force_output();
        }
        break;
      }

      if (males.size() > 1e6) {
        if (verbose) {
          Rcpp::warning("encountered massive population, culling to 1M indiv");
        }
        cul_population(males, 1e6);
      }
      if (females.size() > 1e6) {
        Rcpp::warning("encountered massive population, culling to 1M indiv");

        cul_population(females, 1e6);
      }

      std::array<double, 3> f2 = calc_freq_hawaii< BIRD >(males, females);

      frequencies.add_slice(replicate, t, f2,
                            static_cast<int>( males.size() + females.size() ),
                            males.size(),
                            females.size());
        if (verbose) {
            Rcpp::Rcout << t << " " << f2[0] << "\t" << males.size()  << "\t" <<  females.size() << "\n"; force_output();
        }

    }

    output_pop = males;
    for (const auto& i : females) {
      output_pop.push_back(i);
    }
    return frequencies;
  }

  void additional_death(std::vector< BIRD>& birds,
                        double death_rate,
                        int max_num) {

    if (max_num > birds.size()) max_num = birds.size();

    int num_dead = rndgen.binomial(max_num, death_rate);
    if (num_dead <= 0) return;
    for (int i = 0; i < num_dead; ++i) {
      int index = rndgen.random_number(birds.size());
      birds[index] = birds.back();
      birds.pop_back();
    }
    return;
  }

  void update_pop(std::vector< BIRD >& females,
                  std::vector< BIRD >& males,
                  int number_added,
                  int number_removed) {

    size_t pop_size = females.size() + males.size();

    float frac = females.size() * 1.0 / (pop_size);

    double density_dependent_death_rate = calculate_death_rate(pop_size);

    int males_added = number_added * params.sex_ratio_put;
    int females_added = number_added - males_added;
    if (females_added < 0) females_added = 0;

    update_start_season(females,
                        density_dependent_death_rate,
                        number_removed * frac,
                        females_added);

    update_start_season(males,
                        density_dependent_death_rate,
                        number_removed * (1 - frac),
                        males_added);

    if(females.empty() && males.empty()) {
       return;
    }

    std::vector< BIRD > offspring_male;
    std::vector< BIRD > offspring_female;

    pop_size = males.size() + females.size(); // + offspring_male.size() + offspring_female.size();
    double density_dependent_offspring_rate = calculate_death_rate(pop_size);

    // now, mated females experience additional death
    if (males.size() < females.size()) {
      // not all females are mated, females should be shuffled
      std::shuffle(females.begin(), females.end(), rndgen.rndgen);
    } else {
      // not all males are mated, males should be shuffled
      std::shuffle(males.begin(), males.end(), rndgen.rndgen);
    }

    for (int i = 0; i < females.size() && i < males.size(); ++i) {
      if (rndgen.bernouilli(params.female_death_rate)) {
        females[i] = females.back();
        females.pop_back();
        --i;
      } else {
        generate_offspring(offspring_male,
                           offspring_female,
                           females[i],
                           males[i],
                           density_dependent_offspring_rate,
                           params.clutch_size_mean,
                           params.clutch_size_sd,
                           params.sex_ratio_offspring);
      }
    }

    if (!offspring_male.empty()) {
      males.insert(males.end(),
                   offspring_male.begin(), offspring_male.end());
    }


    if (!offspring_female.empty()) {
      females.insert(females.end(),
                     offspring_female.begin(), offspring_female.end());
    }

    return;
  }

  void generate_offspring(std::vector< BIRD >& offspring_male,
                          std::vector< BIRD >& offspring_female,
                          const BIRD& mama,
                          const BIRD& papa,
                          double offspring_death_rate,
                          int clutch_size,
                          double clutch_sd,
                          double prob_male) {

      if (rndgen.bernouilli(1.0 - params.nest_failure_rate)) {
      // nest is not predated
      int num_offspring = static_cast<int>(rndgen.normal_positive(clutch_size, clutch_sd));

      for (int j = 0; j < num_offspring; ++j) {

        if (rndgen.uniform() > offspring_death_rate) { // immediately check survival to next generation

          BIRD chick(mama.gamete(params.morgan, rndgen),
                     papa.gamete(params.morgan, rndgen),
                     prob_male,
                     rndgen);

          if (chick.get_sex() == female) {
            offspring_female.push_back(std::move(chick));
          } else {
            offspring_male.push_back(std::move(chick));
          }
        }
      }
    }
    return;
  }

  float calculate_death_rate(size_t N) {
    float d = 1.f * N  / params.K;
    float numerator = 1 + powf(1.f * d / params.p, params.b);

    return  1.f - (params.smax + 1.f * (params.smin - params.smax) / (numerator));
  }

  void old_age(std::vector< BIRD >& pop) {
    for (int i = 0; i < pop.size(); ++i) {
      pop[i].age++;
      if (pop[i].age > params.max_age) {
        pop[i] = pop.back();
        pop.pop_back();
        i--;
      }
    }
  }

  void update_start_season(std::vector< BIRD >& input_pop,
                           double death_rate,
                           int number_removed,
                           int number_added) {

    // first, regular death due to old age
    old_age(input_pop);

    // then, death due to survival:
    additional_death(input_pop,
                     death_rate,
                     input_pop.size());

    // then, they can be killed (until number removed is reached)
    if (number_removed > 0) {
      if (number_removed >= input_pop.size()) {
        input_pop.clear();
      } else {
        for (int i = 0; i < number_removed; ++i) {
          int unlucky_indiv = rndgen.random_number(input_pop.size());
          input_pop[unlucky_indiv] = input_pop.back();
          input_pop.pop_back();
          if (input_pop.empty()) break;
        }
      }
    }

    // then we add hawaii individuals
    if (number_added > 0) {
      add_to_population(input_pop,
                        number_added,
                        tag<BIRD>{},
                        input_pop.back().get_sex()
                       );
    }
    return;
  }

  void add_to_population(std::vector<Duck_simple>& population,
                         int number_added, tag<Duck_simple>,
                         const Sex& sex) {
    Duck_simple to_add(1); // Hawaii = 1
    to_add.set_sex(sex);
    for(int i = 0; i < number_added; ++i) {
      population.push_back(to_add);
    }
    return;
  }

  void add_to_population(std::vector<Duck>& population,
                                             int number_added, tag<Duck>,
                                             const Sex& sex) {
    Duck to_add(1); // Hawaii = 1
    to_add.set_sex(sex);
    for(int i = 0; i < number_added; ++i) {
      population.push_back(to_add);
    }
  }

  void add_to_population(std::vector<Duck_emp>& population,
                         int number_added,
                         tag<Duck_emp>,
                         const Sex& sex) {
    for(int i = 0; i < number_added; ++i) {
      Duck_emp to_add = draw_pure();
      to_add.set_sex(sex);
      population.emplace_back(to_add);
    }
  }

  std::vector<Duck_emp> create_base_pop(tag<Duck_emp>) {
    return NumericMatrix_to_emp_bird(base_genomes,
                                     params.empgen,
                                     rndgen);
  }

  std::vector<Duck_simple> create_base_pop(tag<Duck_simple>) {
    return admix();
  }

  std::vector<Duck> create_base_pop(tag<Duck>) {
    BIRD mallard(0.0); // = Mallard;
    BIRD hawaii(1.0); // = Hawaii;

    std::vector< BIRD > population(params.pop_size);

    for(int i = 0; i < params.pop_size; ++i) {
      BIRD parent1 = mallard;
      BIRD parent2 = mallard;

      float freq_hawaii = rndgen.normal_positive(params.frequency_hawaii_duck,
                                                 params.sd_frequency_hawaii);

      if(rndgen.uniform() < freq_hawaii) {
        parent1 = hawaii;
      }

      freq_hawaii = rndgen.normal_positive(params.frequency_hawaii_duck,
                                           params.sd_frequency_hawaii);

      if(rndgen.uniform() < freq_hawaii) {
        parent2 = hawaii;
      }

      double init_prob_female = 0.5;

      population[i] = BIRD(parent1.gamete(params.morgan, rndgen),
                           parent2.gamete(params.morgan ,rndgen),
                           init_prob_female, rndgen);
    }
    return population;
  }

  std::vector< BIRD > admix() {
    BIRD mallard(0.0); // = Mallard;
    BIRD hawaii(1.0); // = Hawaii;

    std::vector< BIRD > population(params.pop_size);

    for(int i = 0; i < params.pop_size; ++i) {
      BIRD parent1 = mallard;
      BIRD parent2 = mallard;

      float freq_hawaii = rndgen.normal_positive(params.frequency_hawaii_duck,
                                                 params.sd_frequency_hawaii);

      if(rndgen.uniform() < freq_hawaii) {
        parent1 = hawaii;
      }

      freq_hawaii = rndgen.normal_positive(params.frequency_hawaii_duck,
                                           params.sd_frequency_hawaii);

      if(rndgen.uniform() < freq_hawaii) {
        parent2 = hawaii;
      }

      double init_prob_female = 0.5;

      population[i] = BIRD(parent1.gamete(params.morgan, rndgen),
                           parent2.gamete(params.morgan ,rndgen),
                           init_prob_female, rndgen);
    }

    for (size_t t = 0; t < params.number_of_generations; ++t) {

      std::vector< BIRD > new_population(params.pop_size);
      for (size_t i = 0; i < params.pop_size; ++i) {
        int index1 = rndgen.random_number(params.pop_size);
        int index2 = rndgen.random_number(params.pop_size);
        while(index2 == index1) index2 = rndgen.random_number(params.pop_size);

        new_population[i] = BIRD(population[index1].gamete(params.morgan, rndgen),
                                 population[index2].gamete(params.morgan, rndgen),
                                 0.5,
                                 rndgen);
      }
      population.swap(new_population);
    }
    return population;
  }

  void cul_population(std::vector<BIRD>& pop, size_t target_size) {
    if (pop.size() < target_size) return;
    pop.erase(pop.begin() + target_size, pop.end());
  }
};


#endif
