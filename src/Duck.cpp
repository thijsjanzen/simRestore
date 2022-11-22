//
//  Duck.cpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Duck.h"
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

Duck::Duck() {
    age = 0;
    freq_hawaii = -1;
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

Duck::Duck(int initLoc)    {
    junction left(0.0, initLoc);
    junction right(1,  -1);
    chromosome1.push_back( left  );
    chromosome1.push_back( right );
    chromosome2.push_back( left  );
    chromosome2.push_back( right );

    freq_hawaii = initLoc;
    age = 0;
    return;
}

Duck::Duck(const std::vector<junction>& c1,
           const std::vector<junction>& c2,
           double prob_female,
           rnd_t& rndgen) :
    chromosome1(c1), chromosome2(c2) {
    calc_freq_hawaii();
    set_nonrandom_sex(prob_female, rndgen);
    age = 0;
}

void Duck::set_nonrandom_sex(double prob_male,
                             rnd_t& rndgen) {
  sex = female;
  if (rndgen.uniform() < prob_male) {
    sex = male;
  }
  return;
}

void Duck::set_random_sex(rnd_t& rndgen) noexcept {
    sex = female;
    if (rndgen.uniform() < 0.5) {
        sex = male;
    }
    return;
}

Duck::Duck(const Duck& other) : age(other.age), chromosome1(other.chromosome1), chromosome2(other.chromosome2),
sex(other.sex), freq_hawaii(other.freq_hawaii)  {
}

Duck& Duck::operator=(const Duck& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
        sex = other.sex;
        freq_hawaii = other.freq_hawaii;
        age = other.age;
    }
    return *this;
}

Duck::Duck(Duck&& other) {
    chromosome1 = other.chromosome1;
    chromosome2 = other.chromosome2;
    sex = other.sex;
    freq_hawaii = other.freq_hawaii;
    age = other.age;
}

Duck& Duck::operator=(Duck&& other) {
    if (this != &other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
        sex = other.sex;
        freq_hawaii = other.freq_hawaii;
        age = other.age;
    }
    return *this;
}

std::vector<junction> Duck::gamete(double morgan, rnd_t& rndgen) const noexcept {
    int num_recom = rndgen.poisson(morgan);
    std::vector<double> recom_pos(num_recom);
    for (int i = 0; i < num_recom; ++i) recom_pos[i] = rndgen.uniform();
    std::sort(recom_pos.begin(), recom_pos.end());
    recom_pos.push_back(1.0);

    std::vector<junction> offspring_chromosome;
    if (rndgen.random_number(2) == 0) {
        offspring_chromosome = recombine(chromosome1, chromosome2, recom_pos);
    } else {
        offspring_chromosome = recombine(chromosome2, chromosome1, recom_pos);
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

void Duck::calc_freq_hawaii() {
    double freq1 = calc_freq_chrom(chromosome1);
    double freq2 = calc_freq_chrom(chromosome2);
    freq_hawaii = 0.5 * (freq1 + freq2);
}

Duck_simple::Duck_simple() {
    freq_hawaii = -1.0;
    age = 0;
}

Duck_simple::Duck_simple(double initLoc) {
    freq_hawaii = initLoc;
    chromosome1 = initLoc;
    chromosome2 = initLoc;
    age = 0;
}

Duck_simple::Duck_simple(double chrom1,
                         double chrom2,
                         double prob_female,
                         rnd_t& rndgen) :
    chromosome1(chrom1), chromosome2(chrom2) {
    freq_hawaii = 0.5 * (chromosome1 + chromosome2);
    set_nonrandom_sex(prob_female, rndgen);
    age = 0;
}

Duck_simple::Duck_simple(const Duck_simple& other)  {
    freq_hawaii = other.get_freq_hawaii();
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    sex = other.sex;
    age = other.age;
}

// DUCK TO DUCK_SIMPLE
Duck_simple::Duck_simple(const Duck& other) {
    sex = other.get_sex();
    chromosome1 = calc_freq_chrom(other.get_chromosome1()); // conversion from std::vector<junction> to double
    chromosome2 = calc_freq_chrom(other.get_chromosome2());
    freq_hawaii = other.get_freq_hawaii();
    age = other.age;
}

Duck_simple& Duck_simple::operator=(const Duck_simple& other) {
    freq_hawaii = other.get_freq_hawaii();
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    sex = other.get_sex();
    age = other.age;
    return *this;
}

void Duck_simple::set_nonrandom_sex(double prob_male,
                                    rnd_t& rndgen) {
  sex = female;
  if (rndgen.uniform() < prob_male) {
    sex = male;
  }
  return;
}

void Duck_simple::set_random_sex(rnd_t& rndgen) noexcept {
    sex = female;
    if (rndgen.uniform() < 0.5) {
        sex = male;
    }
    return;
}

double Duck_simple::gamete(double morgan, rnd_t& rndgen) const noexcept {
    // recombine chromosomes:
    return 0.5 * (chromosome1 + chromosome2);
}


//////////////////////////////////////////////////////////////////////////////
///// DUCK WITH EMPIRICALDATA  ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

Duck_emp::Duck_emp(Duck_emp&& other) {
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    empgen_ = other.get_empgen();
    freq_hawaii = other.get_freq_hawaii();
    sex = other.get_sex();
    age = other.age;
}

Duck_emp& Duck_emp::operator=(Duck_emp&& other) {
    if (this != &other) {
        chromosome1 = other.get_chromosome1();
        chromosome2 = other.get_chromosome2();
        empgen_ = other.get_empgen();
        sex = other.get_sex();
        set_hawaii_freq();
        age = other.age;
    }
    return *this;
}


Duck_emp::Duck_emp(const Duck_emp& other) {
    chromosome1 = other.get_chromosome1();
    chromosome2 = other.get_chromosome2();
    empgen_ = other.get_empgen();
    set_hawaii_freq();
    sex = other.get_sex();
    age = other.age;
}

Duck_emp& Duck_emp::operator=(const Duck_emp& other) {
    if (this != &other) {
        chromosome1 = other.get_chromosome1();
        chromosome2 = other.get_chromosome2();
        empgen_ = other.get_empgen();
        set_hawaii_freq();
        sex = other.get_sex();
        age = other.age;
    }
    return *this;
}

void Duck_emp::set_random_sex(rnd_t& rndgen) {
    if (rndgen.random_number(2)) {
        sex = male;
    } else {
        sex = female;
    }
}

/*
double Duck_emp::calc_anc(const std::vector< double >& allele_freq,
                          double a1,
                          double a2) {

  if (a1 == 0) {
    if((allele_freq[a2 + 4] + allele_freq[a2]) == 0) return 0;

    return(allele_freq[a2 + 4] / (allele_freq[a2 + 4] + allele_freq[a2]));
  }
  if (a2 == 0) {
    if ((allele_freq[a1 + 4] + allele_freq[a1]) == 0) return 0;

    return(allele_freq[a1 + 4] / (allele_freq[a1 + 4] + allele_freq[a1]));
  }

  if (a1 == a2) {
    double denominator = (allele_freq[a1 + 4] + allele_freq[a1]);
    if (denominator == 0) denominator = 1;
    double prob_allele_h = allele_freq[a1 + 4] / denominator;
    return(prob_allele_h * prob_allele_h);
  }
  if (a1 != a2) {
    double denom1 = (allele_freq[a1 + 4] + allele_freq[a1]);
    if (denom1 == 0) denom1 = 1;

    double denom2 = (allele_freq[a2 + 4] + allele_freq[a2]);
    if (denom2 == 0) denom2 = 1;



    double prob_A_h = allele_freq[a1 + 4] / denom1;
    double prob_a_h = allele_freq[a2 + 4] / denom2;;
    double prob_A_m = allele_freq[a1] / denom1;
    double prob_a_m = allele_freq[a2] / denom2;
    return(prob_A_h * prob_a_m + prob_A_m * prob_a_h);
  }

  return 0.0;
                          }*/

double Duck_emp::calc_anc(const std::vector< double >& allele_freq,
                          double a1,
                          double a2) {
  if (a1 == 0) {
    return 0.0;
  }
  if (a2 == 0) {
    return 0.0;
  }

  // prob both alleles mallard:
  double prob_mall = allele_freq[a1] * allele_freq[a2];

 //prob both alleles hawaii:
  double prob_haw  = allele_freq[a1 + 4] * allele_freq[a2 + 4];

  // prob heterozygous
  double prob_het  = allele_freq[a1] * allele_freq[a2 + 4] +
                      allele_freq[a2] * allele_freq[a1 + 4];

  double prob_hawaiian = prob_het * 0.5 / (prob_mall + prob_het + prob_haw) +
                         prob_haw / (prob_mall + prob_het + prob_haw);

//  Rcpp::Rcout << allele_freq[0] << " " << a1 << " " << a2 << " " << prob_mall << " " << prob_haw << " " << prob_het << " " << prob_hawaiian << "\n";
//  R_FlushConsole();
//  R_ProcessEvents();
//  R_CheckUserInterrupt();


  return(prob_hawaiian);
}


double Duck_emp::set_hawaii_freq() {
  //  double freq_hawaii = 0.0;
  std::vector<double> freqs(empgen_.anc_info.size(), 0.0);
    for (size_t i = 0; i < empgen_.anc_info.size(); ++i) {
        int pos = empgen_.anc_info[i][0] - 1; // R indexing!
        freqs[i] = calc_anc(empgen_.anc_info[i],
                                chromosome1[pos],
                                chromosome2[pos]);
    }
  double freq_hawaii = std::accumulate(freqs.begin(), freqs.end(), 0.0);
  double output = freq_hawaii * 1.0 / (empgen_.anc_info.size());
  return output;
}

std::tuple< std::vector<int>, emp_genome> Duck_emp::gamete(double morgan,
                                                           rnd_t& rndgen) const {

    std::vector<size_t> recom_pos = empgen_.recompos(morgan,
                                                     rndgen);

    if (recom_pos.size() == 1) {
        if(rndgen.random_number(2)) {
            return std::make_tuple(chromosome1, empgen_);
        }
        return std::make_tuple(chromosome2, empgen_);
    }

    std::vector < std::vector<int>::const_iterator > iters = {chromosome1.begin(),
                                                              chromosome2.begin()};
    std::vector< int > recombined_chromosome;
    int index = rndgen.random_number(2);
    size_t prev_start = 0;

    for(size_t i = 0; i < recom_pos.size(); ++i) {
        auto start = iters[index] + prev_start;
        auto end   = iters[index] + recom_pos[i];

        prev_start = recom_pos[i];
        recombined_chromosome.insert(recombined_chromosome.end(), start, end);
        index = 1 - index;
    }

    return std::make_tuple(recombined_chromosome, empgen_);
}
