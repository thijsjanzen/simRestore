#pragma once

#include "rand_t.h"
#include <tuple>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>

struct junction {
    long double pos;
    int right;

    junction() {};
    junction(long double loc, int B);
    junction(const junction& other);
    junction& operator=(const junction& other);
};

enum Sex {male, female};

using chrom = std::vector< junction >;
using genome = std::vector< chrom >;

struct organism {
    organism();
    organism(double init_freq, size_t num_chromosomes);
    organism(const genome& c1,
             const genome& c2,
             double prob_male,
             rnd_t& rndgen);


    organism(const organism& other);
    organism& operator=(const organism& other);

    organism(organism&& other);
    organism& operator=(organism&& other);

    void set_nonrandom_sex(double prob_male, rnd_t& rndgen);
    void set_sex(Sex s) {sex = s;}

    genome gamete(const std::vector<double>& morgan,
                  rnd_t& rndgen) const noexcept;

    const genome& get_chromosome1() const noexcept {return chromosome1;}
    const genome& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_anc() const noexcept {return freq_anc;}
    const Sex& get_sex() const noexcept {return sex;}
    int age;

 private:
    genome chromosome1;
    genome chromosome2;
    Sex sex;
    double freq_anc;
    void calc_freq_anc();
};

struct organism_simple {
    organism_simple();
    organism_simple(double initLoc, size_t num_chromosomes);

    organism_simple(double chrom1, double chrom2,
                    double prob_male, rnd_t& rndgen);
    organism_simple(const organism_simple& other);
    organism_simple& operator=(const organism_simple& other);

    void set_nonrandom_sex(double prob_male, rnd_t& rndgen);
    void set_sex(Sex s) {sex = s;}

    double gamete(const std::vector<double>& morgan,
                  rnd_t& rndgen) const noexcept;

    const double& get_chromosome1() const noexcept {return chromosome1;}
    const double& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_anc() const noexcept {return freq_anc;}
    const Sex& get_sex() const noexcept {return sex;}

    int age;

 private:
    double chromosome1;
    double chromosome2;
    double freq_anc;
    Sex sex;
};

struct emp_genome {
  // placeholder
};

struct organism_emp {
  // placeholder;
  void set_sex(bool s);
};
