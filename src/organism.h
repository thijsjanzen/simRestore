#pragma once

#include <tuple>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "rand_t.h"


struct junction {
    long double pos;
    int right;

    junction();
    junction(long double loc, int B) ;
    junction(const junction& other);
    bool operator <(const junction& other) const noexcept;
    junction& operator=(const junction& other);
};

enum Sex {male, female};

struct organism {

    organism();
    organism(int initLoc);
    organism(const std::vector<junction>& c1,
         const std::vector<junction>& c2,
         double prob_male,
         rnd_t& rndgen);


    organism(const organism& other);
    organism& operator=(const organism& other);

    organism(organism&& other);
    organism& operator=(organism&& other);

    void set_nonrandom_sex(double prob_male, rnd_t& rndgen);
    void set_random_sex(rnd_t& rdngen) noexcept;
    void set_sex(Sex s) {sex = s;}

    std::vector<junction> gamete(double morgan, rnd_t& rndgen) const noexcept;

    const std::vector< junction >& get_chromosome1() const noexcept {return chromosome1;}
    const std::vector< junction >& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_hawaii() const noexcept {return freq_hawaii;}
    const Sex& get_sex() const noexcept {return sex;}
    int age;

private:
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;
    Sex sex;
    double freq_hawaii;
    void calc_freq_hawaii();
};

struct organism_simple {

    organism_simple();
    organism_simple(double initLoc);

    organism_simple(const organism& other);

    organism_simple(double chrom1, double chrom2, double prob_male, rnd_t& rndgen);
    organism_simple(const organism_simple& other);
    organism_simple& operator=(const organism_simple& other);

    void set_nonrandom_sex(double prob_male, rnd_t& rndgen);
    void set_random_sex(rnd_t& rndgen) noexcept;
    void set_sex(Sex s) {sex = s;}

    double gamete(double morgan, rnd_t& rndgen) const noexcept;

    const double& get_chromosome1() const noexcept {return chromosome1;}
    const double& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_hawaii() const noexcept {return freq_hawaii;}
    const Sex& get_sex() const noexcept {return sex;}

    int age;
private:
    double chromosome1;
    double chromosome2;
    double freq_hawaii;
    Sex sex;
};

struct emp_genome {
  // placeholder
};

struct organism_emp {
  // placeholder;
  void set_sex(bool s);
};
