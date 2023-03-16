#ifndef MAIN_HPP
#define MAIN_HPP

#include "organism.h"
#include <vector>
#include <array>
#include <Rcpp.h>
using namespace Rcpp;


struct output_entry {
  int replicate;
  size_t t;
  double frequency;
  double frequency_males;
  double frequency_females;
  size_t pop_size;
  size_t num_males;
  size_t num_females;

  output_entry() {

  }

  output_entry(int r, size_t time, double f, double f_m, double f_f,
               size_t n, size_t n_m, size_t n_f) :
  replicate(r),
  t(time),
  frequency(f),
  frequency_males(f_m),
  frequency_females(f_f),
  pop_size(n),
  num_males(n_m),
  num_females(n_f) {
  }

  output_entry(const output_entry& other) {
    replicate = other.replicate;
    t = other.t;
    frequency = other.frequency;
    frequency_males = other.frequency_males;
    frequency_females = other.frequency_females;
    pop_size = other.pop_size;
    num_males = other.num_males;
    num_females = other.num_females;
  }
};


struct output_data {

  output_data() {}

  void add_slice(int replicate,
                 size_t t,
                 std::array<double, 3> f,
                 size_t n,
                 size_t n_m,
                 size_t n_f) {
    data_.emplace_back(output_entry(replicate, t, f[0], f[1], f[2], n, n_m, n_f));
  }

  const std::vector< output_entry>& data() const noexcept {return data_;}
  const int size() const noexcept {return data_.size();}

  const output_entry operator[](size_t i){ return data_[i];}

private:
  std::vector< output_entry > data_;
};



void force_output();

#endif
