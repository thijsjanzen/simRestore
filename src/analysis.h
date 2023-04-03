#pragma once

#include <vector>
#include <algorithm>

#include "rand_t.h"
#include "main_functions.h"

template <typename> struct tag { };

struct parameters {

  parameters() {
    pop_size = 100;
    starting_freq = 0.2;
    sd_starting_freq = 0.05;
    number_of_generations = 20;
    K = 400;
    morgan = {1.0};
    female_death_rate = 0.2;
    male_death_rate = 0.0;
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
    sex_ratio_pull = 0.5;
    sex_ratio_offspring = 0.5;
    put_ancestry = 1.0;
  }

  parameters(int p,
             float f,
             float s_f,
             int num_gen,
             int k,
             std::vector<double> m,
             double f_d_r,
             double m_d_r,
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
             double sex_ratio_add,
             double sex_ratio_remove,
             double sex_ratio_o,
             double put_anc) :
  pop_size(p),
  starting_freq(f),
  sd_starting_freq(s_f),
  number_of_generations(num_gen),
  K(k),
  morgan(m),
  female_death_rate(f_d_r),
  male_death_rate(m_d_r),
  nest_failure_rate(n_f_r),
  establishment_burnin(burnin),
  max_age(m_a),
  clutch_size_mean(cs_m),
  clutch_size_sd(cs_sd),
  smin(smin_in),
  smax(smax_in),
  p(p_in),
  b(b_in),
  sex_ratio_put(sex_ratio_add),
  sex_ratio_pull(sex_ratio_remove),
  sex_ratio_offspring(sex_ratio_o),
  put_ancestry(put_anc),
  empgen(emp_gen_in){
  }

  int pop_size;
  float starting_freq;
  float sd_starting_freq;
  int number_of_generations;
  int K;
  std::vector<double> morgan;
  double female_death_rate;
  double male_death_rate;
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
  double sex_ratio_pull;
  double sex_ratio_offspring;

  double put_ancestry; // ancestry used when putting

  emp_genome empgen; // only used in molecular analysis
};

template <class ANIMAL>
std::array<double, 3> calc_freq_focal(std::vector< ANIMAL >& f,
                                      std::vector< ANIMAL >& m) {
  std::array<double, 3> avg_freq = {0.0, 0.0, 0.0};

  for (auto& i : f) {
    double freq_i = i.get_freq_anc();
    avg_freq[2] += freq_i;
    avg_freq[0] += freq_i;
  }
  for (auto& i : m) {
    double freq_i = i.get_freq_anc();
    avg_freq[1] += freq_i;
    avg_freq[0] += freq_i;
  }

  avg_freq[0] *= 1.0 / (f.size() + m.size());
  avg_freq[1] *= 1.0 / (m.size());
  avg_freq[2] *= 1.0 / (f.size());

  return avg_freq;
}


template <typename ANIMAL>
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
  std::vector< organism_emp > pure;
  std::vector< ANIMAL > output_pop;

  std::vector<int> recorded_death_ages;

  // member functions
  void increase_replicate() {
    replicate++;
  }

  void set_pure(const std::vector<organism_emp>& p) {
    pure = p;
  }

  organism_emp draw_pure() {
    int index = rndgen.random_number(static_cast<int>(pure.size()));
    return pure[index];
  }

  output_data do_analysis() {
    std::vector< ANIMAL > base_pop = create_base_pop(tag<ANIMAL>{});

    output_data run_result = simulate_policy(base_pop);

    return run_result;
  }

private:
  int replicate;

  output_data simulate_policy(const std::vector< ANIMAL >& base_pop) {
    std::vector< ANIMAL > males;
    std::vector< ANIMAL > females;

    int num_females = static_cast<int>((base_pop.size() / 2.0) * 1.0);

    for(int i = 0; i < num_females; ++i) {
      females.push_back(base_pop[i]);
    }
    for(int i = num_females; i < base_pop.size(); ++i) {
      males.push_back(base_pop[i]);
    }

    output_data frequencies;
    std::array<double, 3> f1 = calc_freq_focal<ANIMAL>(females, males);

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

      std::array<double, 3> f2 = calc_freq_focal< ANIMAL >(males, females);

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

  void additional_death(std::vector< ANIMAL>& ANIMALs,
                        double death_rate,
                        int max_num) {

    if (max_num > ANIMALs.size()) max_num = ANIMALs.size();

    int num_dead = rndgen.binomial(max_num, death_rate);
    if (num_dead <= 0) return;
    for (int i = 0; i < num_dead; ++i) {
      int index = rndgen.random_number(ANIMALs.size());
      recorded_death_ages.push_back(ANIMALs[index].age);
      ANIMALs[index] = ANIMALs.back();
      ANIMALs.pop_back();
    }
    return;
  }

  void update_pop(std::vector< ANIMAL >& females,
                  std::vector< ANIMAL >& males,
                  int number_added,
                  int number_removed) {

    size_t pop_size = females.size() + males.size();

    double density_dependent_death_rate = calculate_death_rate(pop_size);

    int males_added = number_added * params.sex_ratio_put;
    int females_added = number_added - males_added;
    if (females_added < 0) females_added = 0;

    update_start_season(females,
                        density_dependent_death_rate,
                        number_removed * (1.0 - params.sex_ratio_pull),
                        females_added);

    update_start_season(males,
                        density_dependent_death_rate,
                        number_removed * params.sex_ratio_pull,
                        males_added);

    if(females.empty() && males.empty()) {
       return;
    }

    std::vector< ANIMAL > offspring_male;
    std::vector< ANIMAL > offspring_female;

    pop_size = males.size() + females.size();
    double density_dependent_offspring_rate = calculate_death_rate(pop_size);


    if (males.size() < females.size()) {
      // not all females are mated, females should be shuffled
      std::shuffle(females.begin(), females.end(), rndgen.rndgen);
    } else {
      // not all males are mated, males should be shuffled
      std::shuffle(males.begin(), males.end(), rndgen.rndgen);
    }

    for (int i = 0, j = 0; i < females.size() && j < males.size(); ++i, ++j) {
      // now, mated females and females experience additional death
      if (rndgen.bernouilli(params.female_death_rate)) {
        recorded_death_ages.push_back(females[i].age);
        females[i] = females.back();
        females.pop_back();
        --i;
      } else {
        generate_offspring(offspring_male,
                           offspring_female,
                           females[i],
                           males[j],
                           density_dependent_offspring_rate,
                           params.clutch_size_mean,
                           params.clutch_size_sd,
                           params.sex_ratio_offspring);

        if (rndgen.bernouilli(params.male_death_rate)) {
          recorded_death_ages.push_back(males[j].age);
          males[j] = males.back();
          males.pop_back();
          --j;
        }
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

  void generate_offspring(std::vector< ANIMAL >& offspring_male,
                          std::vector< ANIMAL >& offspring_female,
                          const ANIMAL& mama,
                          const ANIMAL& papa,
                          double offspring_death_rate,
                          int clutch_size,
                          double clutch_sd,
                          double prob_male) {

      if (rndgen.bernouilli(1.0 - params.nest_failure_rate)) {
      // nest is not predated
      int num_offspring = static_cast<int>(rndgen.normal_positive(clutch_size, clutch_sd));

      for (int j = 0; j < num_offspring; ++j) {

        if (rndgen.uniform() > offspring_death_rate) { // immediately check survival to next generation

          ANIMAL chick(mama.gamete(params.morgan, rndgen),
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

  void old_age(std::vector< ANIMAL >& pop) {
    for (int i = 0; i < pop.size(); ++i) {
      pop[i].age++;
      if (pop[i].age > params.max_age) {
        recorded_death_ages.push_back(pop[i].age);
        pop[i] = pop.back();
        pop.pop_back();
        i--;
      }
    }
  }

  void update_start_season(std::vector< ANIMAL >& input_pop,
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

    // then we add pure individuals
    if (number_added > 0) {
      add_to_population(input_pop,
                        number_added,
                        tag<ANIMAL>{},
                        input_pop.back().get_sex()
                       );
    }
    return;
  }

  void add_to_population(std::vector<organism_simple>& population,
                         int number_added, tag<organism_simple>,
                         const Sex& sex) {
    organism_simple to_add(params.put_ancestry, params.morgan.size());
    to_add.set_sex(sex);
    for(int i = 0; i < number_added; ++i) {
      population.push_back(to_add);
    }
    return;
  }

  void add_to_population(std::vector<organism>& population,
                         int number_added, tag<organism>,
                         const Sex& sex) {
    organism to_add(params.put_ancestry, params.morgan.size());
    to_add.set_sex(sex);
    for(int i = 0; i < number_added; ++i) {
      population.push_back(to_add);
    }
  }

  void add_to_population(std::vector<organism_emp>& population,
                         int number_added,
                         tag<organism_emp>,
                         const Sex& sex) {
    for(int i = 0; i < number_added; ++i) {
      organism_emp to_add = draw_pure();
      to_add.set_sex(sex);
      population.emplace_back(to_add);
    }
  }

  std::vector<organism_simple> create_base_pop(tag<organism_simple>) {
    return admix();
  }

  std::vector<organism> create_base_pop(tag<organism>) {
    ANIMAL base_indiv(0.0, params.morgan.size());
    ANIMAL target_indiv(1.0, params.morgan.size());

    std::vector< ANIMAL > population(params.pop_size);

    for(int i = 0; i < params.pop_size; ++i) {
      ANIMAL parent1 = base_indiv;
      ANIMAL parent2 = base_indiv;

      float freq_focal = rndgen.normal_positive(params.starting_freq,
                                                 params.sd_starting_freq);

      if(rndgen.uniform() < freq_focal) {
        parent1 = target_indiv;
      }

      freq_focal = rndgen.normal_positive(params.starting_freq,
                                           params.sd_starting_freq);

      if(rndgen.uniform() < freq_focal) {
        parent2 = target_indiv;
      }

      double init_prob_female = 0.5;

      population[i] = ANIMAL(parent1.gamete(params.morgan, rndgen),
                           parent2.gamete(params.morgan ,rndgen),
                           init_prob_female, rndgen);
    }
    return population;
  }

  std::vector< ANIMAL > admix() {
    ANIMAL base_indiv(0.0, params.morgan.size());
    ANIMAL target_indiv(1.0, params.morgan.size());

    std::vector< ANIMAL > population(params.pop_size);

    for(int i = 0; i < params.pop_size; ++i) {
      ANIMAL parent1 = base_indiv;
      ANIMAL parent2 = base_indiv;

      float freq_focal = rndgen.normal_positive(params.starting_freq,
                                                 params.sd_starting_freq);

      if(rndgen.uniform() < freq_focal) {
        parent1 = target_indiv;
      }

      freq_focal = rndgen.normal_positive(params.starting_freq,
                                           params.sd_starting_freq);

      if(rndgen.uniform() < freq_focal) {
        parent2 = target_indiv;
      }

      double init_prob_female = 0.5;

      population[i] = ANIMAL(parent1.gamete(params.morgan, rndgen),
                           parent2.gamete(params.morgan ,rndgen),
                           init_prob_female, rndgen);
    }

    for (size_t t = 0; t < params.number_of_generations; ++t) {

      std::vector< ANIMAL > new_population(params.pop_size);
      for (size_t i = 0; i < params.pop_size; ++i) {
        int index1 = rndgen.random_number(params.pop_size);
        int index2 = rndgen.random_number(params.pop_size);
        while(index2 == index1) index2 = rndgen.random_number(params.pop_size);

        new_population[i] = ANIMAL(population[index1].gamete(params.morgan, rndgen),
                                 population[index2].gamete(params.morgan, rndgen),
                                 0.5,
                                 rndgen);
      }
      population.swap(new_population);
    }
    return population;
  }

  void cul_population(std::vector<ANIMAL>& pop, size_t target_size) {
    if (pop.size() < target_size) return;
    pop.erase(pop.begin() + target_size, pop.end());
  }
};

