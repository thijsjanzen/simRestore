// Copyright 2023 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#pragma once

#include <stdio.h>
#include <utility>
#include <algorithm>
#include <tuple>
#include <vector>
#include "rand_t.h" // NOLINT [build/include_subdir]

struct junction {
    long double pos;
    int right;

    junction() {}
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
             rnd_t* rndgen);


    organism(const organism& other);
    organism& operator=(const organism& other);

    organism(organism&& other);
    organism& operator=(organism&& other);

    void set_nonrandom_sex(double prob_male, rnd_t* rndgen);
    void set_sex(Sex s) {sex = s;}

    genome gamete(const std::vector<double>& morgan,
                  rnd_t* rndgen) const noexcept;

    const genome& get_chromosome1() const noexcept {return chromosome1;}
    const genome& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_anc() const noexcept {return freq_anc;}
    const Sex& get_sex() const noexcept {return sex;}
    int age;

    std::vector<std::vector<double>> get_genomic_info(int t, int replicate,
                                                      int indiv) const;

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

    organism_simple(const std::vector<double>& chrom1,
                    const std::vector<double>& chrom2,
                    double prob_male, rnd_t* rndgen);

    organism_simple(const organism_simple& other);
    organism_simple& operator=(const organism_simple& other);

    void set_nonrandom_sex(double prob_male, rnd_t* rndgen);
    void set_sex(Sex s) {sex = s;}

    std::vector<double>  gamete(const std::vector<double>& morgan,
                  rnd_t* rndgen) const noexcept;

    const std::vector<double>& get_chromosome1() const noexcept {return chromosome1;}
    const std::vector<double>& get_chromosome2() const noexcept {return chromosome2;}
    const double& get_freq_anc() const noexcept {return freq_anc;}
    const Sex& get_sex() const noexcept {return sex;}

    int age;

    std::vector<std::vector<double>> get_genomic_info(int t, int replicate,
                                                      int indiv) const;

 private:
    std::vector<double> chromosome1;
    std::vector<double> chromosome2;
    double freq_anc;
    Sex sex;

    double calc_freq_chrom(const std::vector<double>& chrom);
};

struct emp_genome {
  // placeholder
};

struct organism_emp {
  // placeholder;
  void set_sex(bool s);
};
