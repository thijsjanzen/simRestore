//
//  organism.cpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "organism.h"
#include <algorithm>
#include <cassert>

#include "Rcpp.h"

std::vector<junction> recombine(const std::vector<junction>& chromosome1,
                                const std::vector<junction>& chromosome2,
                                const std::vector<double>& recom_positions) {

    static auto tl_go = decltype(chromosome1){};
    assert(!chromosome1.empty());    // not strictly enforced by code
    assert(!chromosome2.empty());    // not strictly enforced bu code

    // we need something that is cheaply swappable:
    auto* g1 = &chromosome1;
    auto* g2 = &chromosome2;
    auto& go = tl_go;   // offspring genome: recycle what's already there...
    go.clear();

    // predicate for lower_bound
    auto less = [](const auto& j, double p) { return j.pos < p; };

    // helper lambda to get the value just *before* it.
    // we store the value to the right of a recombination-point but we need the value to the left:
    auto value_at = [](auto begin, auto it) { return (begin != it) ? (it - 1)->right : -1; };

    double left_pos = 0.0;
    auto go_val = -1;
    for (auto right_pos : recom_positions) {
        auto it = std::lower_bound(g1->cbegin(), g1->cend(), left_pos, less);
        auto last = std::lower_bound(it, g1->cend(), right_pos, less);
        // [g1.first, it) : part of the genome *before* left_pos.
        // [it, last) : part of the genome *after or equal to* left_pos but *before* right_pos.
        auto g1_val = value_at(g1->cbegin(), it);
        if (g1_val != go_val) {
            if (it == last || it->pos != left_pos) {
                go.emplace_back(left_pos, g1_val);   // insert change to match
            }
            else {
                ++it;    // corner case: skip spurious double-change
            }
        }
        go.insert(go.end(), it, last);      // append [it, last)
        go_val = value_at(go.begin(), go.end());
        std::swap(g1, g2);
        left_pos = right_pos;
    }
    go.emplace_back(junction(1.0, -1));
    return go;
}

organism::organism() {
    age = 0;
    freq_anc = -1;
}

junction::junction() {

}

junction::junction(long double loc, int B)  {
    pos = loc;
    right = B;
}

junction::junction(const junction& other) {
    pos = other.pos;
    right = other.right;
}

bool junction::operator <(const junction& other) const noexcept {
    return(pos < other.pos);
}

junction& junction::operator=(const junction& other) {
    pos = other.pos;
    right = other.right;
    return *this;
}

organism::organism(double init_freq, size_t num_chromosomes)    {
  for (size_t i = 0; i < num_chromosomes; ++i) {

      chrom chrom1;
      chrom chrom2;

      if (init_freq == 1.0) {
        junction left(0.0, 1);
        junction right(1,  -1);
        chrom1.push_back( left  );
        chrom1.push_back( right );
        chrom2.push_back( left  );
        chrom2.push_back( right );

        freq_anc = init_freq;
        age = 0;
      } else if (init_freq == 0.0) {
        junction left(0.0, 0);
        junction right(1,  -1);
        chrom1.push_back( left  );
        chrom1.push_back( right );
        chrom2.push_back( left  );
        chrom2.push_back( right );

        freq_anc = init_freq;
        age = 0;
      } else {
        // create two alternative chromosomes:
        junction left(0.0, 0);
        junction middle(init_freq, 1.0);
        junction right(1,  -1);

        junction left2(0.0, 1);
        junction middle2(1.0 - init_freq, 0.0);
        junction right2(1,  -1);

        chrom1 = {left, middle, right};
        chrom2 = {left2, middle2, right2};
      }
      chromosome1.push_back(chrom1);
      chromosome2.push_back(chrom2);
    }

    freq_anc = init_freq;
    age = 0;
    return;
}

organism::organism(const genome& c1,
                   const genome& c2,
                   double prob_female,
                   rnd_t& rndgen) :
    chromosome1(c1), chromosome2(c2) {
    calc_freq_anc();
    set_nonrandom_sex(prob_female, rndgen);
    age = 0;
}

void organism::set_nonrandom_sex(double prob_male,
                             rnd_t& rndgen) {
  sex = female;
  if (rndgen.uniform() < prob_male) {
    sex = male;
  }
  return;
}

void organism::set_random_sex(rnd_t& rndgen) noexcept {
    sex = female;
    if (rndgen.uniform() < 0.5) {
        sex = male;
    }
    return;
}

organism::organism(const organism& other) : age(other.age), chromosome1(other.chromosome1), chromosome2(other.chromosome2),
sex(other.sex), freq_anc(other.freq_anc)  {
}

organism& organism::operator=(const organism& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
        sex = other.sex;
        freq_anc = other.freq_anc;
        age = other.age;
    }
    return *this;
}

organism::organism(organism&& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
    sex = other.sex;
    freq_anc = other.freq_anc;
    age = other.age;
}

organism& organism::operator=(organism&& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
        sex = other.sex;
        freq_anc = other.freq_anc;
        age = other.age;
    }
    return *this;
}

genome organism::gamete(const std::vector<double>& morgan,
                        rnd_t& rndgen) const noexcept {

    genome offspring_chromosome;
    for (size_t m = 0; m < morgan.size(); ++m) {
      int num_recom = rndgen.poisson(morgan[m]);
      std::vector<double> recom_pos(num_recom);
      for (int i = 0; i < num_recom; ++i) recom_pos[i] = rndgen.uniform();
      std::sort(recom_pos.begin(), recom_pos.end());
      recom_pos.push_back(1.0);

      chrom offspring_chrom;
      if (rndgen.random_number(2) == 0) {
        offspring_chrom = recombine(chromosome1[m], chromosome2[m], recom_pos);
      } else {
        offspring_chrom = recombine(chromosome2[m], chromosome1[m], recom_pos);
      }
      offspring_chromosome.push_back(offspring_chrom);
    }

    return offspring_chromosome;
}

double calc_freq_chrom(const std::vector< junction >& chrom) {
    double freq = 0.0;
    if (chrom.size() < 2) return 0.0;

    for (size_t i = 1; i < chrom.size(); ++i) {
        double stretch = chrom[i].pos - chrom[i - 1].pos;
        freq += stretch * chrom[i - 1].right;
    }
    return freq;
}

double calc_freq_genome(const genome& g) {
  double freq = 0.0;
  for (const auto& i : g) {
    freq += calc_freq_chrom(i);
  }
  freq *= 1.0 / g.size();
  return freq;

}

void organism::calc_freq_anc() {
    double freq1 = calc_freq_genome(chromosome1);
    double freq2 = calc_freq_genome(chromosome2);
    freq_anc = 0.5 * (freq1 + freq2);
}

organism_simple::organism_simple() {
    freq_anc = -1.0;
    age = 0;
}

organism_simple::organism_simple(double initLoc, size_t num_chromosomes) {
    freq_anc = initLoc;
    chromosome1 = initLoc;
    chromosome2 = initLoc;
    age = 0;
}

organism_simple::organism_simple(double chrom1,
                         double chrom2,
                         double prob_female,
                         rnd_t& rndgen) :
    chromosome1(chrom1), chromosome2(chrom2) {
    freq_anc = 0.5 * (chromosome1 + chromosome2);
    set_nonrandom_sex(prob_female, rndgen);
    age = 0;
}

organism_simple::organism_simple(const organism_simple& other)  {
    freq_anc = other.get_freq_anc();
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    sex = other.sex;
    age = other.age;
}

// organism TO organism_SIMPLE
organism_simple::organism_simple(const organism& other) {
    sex = other.get_sex();
    chromosome1 = calc_freq_genome(other.get_chromosome1()); // conversion from genome to double
    chromosome2 = calc_freq_genome(other.get_chromosome2());
    freq_anc = other.get_freq_anc();
    age = other.age;
}

organism_simple& organism_simple::operator=(const organism_simple& other) {
    freq_anc = other.get_freq_anc();
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    sex = other.get_sex();
    age = other.age;
    return *this;
}

void organism_simple::set_nonrandom_sex(double prob_male,
                                    rnd_t& rndgen) {
  sex = female;
  if (rndgen.uniform() < prob_male) {
    sex = male;
  }
  return;
}

void organism_simple::set_random_sex(rnd_t& rndgen) noexcept {
    sex = female;
    if (rndgen.uniform() < 0.5) {
        sex = male;
    }
    return;
}

double organism_simple::gamete(const std::vector<double>& morgan, rnd_t& rndgen) const noexcept {
    // recombine chromosomes:
    return 0.5 * (chromosome1 + chromosome2);
}

