//
//  Fish.hpp
//
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef Duck_hpp
#define Duck_hpp

#include <tuple>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "rand_t.h"

struct emp_genome {
    std::vector< double > cdf_;
    std::vector< std::vector< double > > anc_info;

    emp_genome() {
    }

    emp_genome(const emp_genome& other) {
        cdf_ = other.cdf_;
        anc_info = other.anc_info;
    }

    emp_genome& operator=(const emp_genome& other) {
        if (this != &other) {
            cdf_ = other.cdf_;
            anc_info = other.anc_info;
        }
        return *this;
    }

    template <typename T>
    emp_genome(const std::vector<T>& positions,
               const std::vector< std::vector< double >>& anc) : anc_info(anc) {
        double total_sum = std::accumulate(positions.begin(),
                                           positions.end(), 0.0);
        double s = 0.0;
        double mult = 1.0 / total_sum;
        cdf_.resize(positions.size());
        for (size_t i = 0; i < positions.size(); ++i) {
            s += positions[i] * mult;
            cdf_[i] = s;
        }

        return;
    }

    emp_genome(const std::vector< std::vector< double >>& anc) : anc_info(anc)
    {}

    emp_genome(const std::vector<double>& positions) {
        double total_sum = std::accumulate(positions.begin(),
                                           positions.end(), 0.0);
        double s = 0.0;
        double mult = 1.0 / total_sum;
        cdf_.resize(positions.size());
        for (size_t i = 0; i < positions.size(); ++i) {
            s += positions[i] * mult;
            cdf_[i] = s;
        }
        return;
    }

    size_t index_from_cdf(double p) const {
        // find index belonging to p
        return static_cast<size_t>(std::distance(cdf_.begin(),
                                                 std::lower_bound(cdf_.begin(),
                                                                  cdf_.end(),
                                                                  p)));
    }

    std::vector< size_t > recompos(double morgan,
                                   rnd_t& rndgen) const {
        size_t num_break_points = rndgen.poisson(morgan);
        std::vector< size_t > indices;
        for(size_t i = 0; i < num_break_points; ++i) {
            auto found_index = index_from_cdf(rndgen.uniform());
            if (found_index > 0) {
                indices.push_back(found_index);
            }
        }
        std::sort(indices.begin(), indices.end());
        indices.push_back(cdf_.size());
        return indices;
    }
};

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

struct Duck {

    Duck();
    Duck(int initLoc);
    Duck(const std::vector<junction>& c1,
         const std::vector<junction>& c2,
         double prob_male,
         rnd_t& rndgen);


    Duck(const Duck& other);
    Duck& operator=(const Duck& other);

    Duck(Duck&& other);
    Duck& operator=(Duck&& other);

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

struct Duck_simple {

    Duck_simple();
    Duck_simple(double initLoc);

    Duck_simple(const Duck& other);

    Duck_simple(double chrom1, double chrom2, double prob_male, rnd_t& rndgen);
    Duck_simple(const Duck_simple& other);
    Duck_simple& operator=(const Duck_simple& other);

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

struct Duck_emp {

    Duck_emp()
    {freq_hawaii = -1.0; age = 0;}

    Duck_emp(const std::tuple< std::vector< int >, emp_genome >& p1,
             const std::tuple< std::vector< int >, emp_genome >& p2,
             rnd_t& rndgen) {
        chromosome1 = std::get<0>(p1);
        chromosome2 = std::get<0>(p2);
        empgen_ = std::get<1>(p1);
        set_random_sex(rndgen);
        set_hawaii_freq();
        age = 0;
    }

    Duck_emp(Duck_emp&& other);
    Duck_emp& operator=(Duck_emp&& other);
    Duck_emp(const Duck_emp& other);
    Duck_emp& operator=(const Duck_emp& other);

    std::tuple< std::vector<int>, emp_genome> gamete(double morgan,
                                                     rnd_t& rndgen) const;

    const std::vector< int >& get_chromosome1() const noexcept {return chromosome1;}
    const std::vector< int >& get_chromosome2() const noexcept {return chromosome2;}
    const emp_genome& get_empgen() const noexcept {return empgen_;}
    const Sex& get_sex() const noexcept {return sex;}
    double get_freq_hawaii()  {
        //set_hawaii_freq();
        //return freq_hawaii;.
      return set_hawaii_freq();
    }

    void set_random_sex(rnd_t& rndgen);
    void set_sex(Sex s) {sex = s;}
    int age;
private:
    std::vector< int > chromosome1;
    std::vector< int > chromosome2;
    Sex sex;
    double freq_hawaii;
    emp_genome empgen_;

    double set_hawaii_freq();
    double calc_anc(const std::vector< double >& allele_freq,
                    double allele1,
                    double allele2);
};


#endif /* Duck_hpp */
